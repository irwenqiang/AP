#include <string>
#include <vector>
#include <iostream>
#include <limits>
#include <algorithm>
#include <graphlab.hpp>
#include <boost/spirit/include/qi.hpp>  
#include <boost/spirit/include/phoenix_core.hpp>    
#include <boost/spirit/include/phoenix_operator.hpp> 
#include <boost/spirit/include/phoenix_stl.hpp>
#include <boost/unordered_set.hpp>

using namespace std;

float former_maxk;
bool perform_scatter;

/*
 *This is the message exchange between vertexs
 */
struct center_msg : public graphlab::IS_POD_TYPE {

	int target_vertex;
	//0: max(a(i, k') + s(i, k')),	s.t. k'!= k
	float max_value;
	//min(0, r(k, k) + sum(max(0, r(i', k))))	s.t. i'!= i
	float min_value;

	center_msg() {}
	explicit center_msg(int tv, float max, float min): target_vertex(tv), max_value(max), min_value(min) {}
	
};

/*
 *Vertex_data_type:
 *vertex_id:int
 *similirity s:float
 */
struct vertex_data : public graphlab::IS_POD_TYPE {

	int vertex_id;
	std::vector<center_msg> msg;
	
	float s;
	float r;
	float a;

	vertex_data() { }
	explicit vertex_data(int id, float fs, float fr, float fa):vertex_id(id), s(fs), r(fr), a(fa) { }

};

/*
 *Edge_data_type:
 *similirity s:vector<float>
 *responsibility r:vector<float>
 *avaliability   a:vector<float>
 */

struct edge_data : public graphlab::IS_POD_TYPE {

	float s;
	float r;
	float a;

	edge_data(){ }
	explicit edge_data(float vs, float vr, float va):s(vs), r(vr), a(va) { }

};

typedef graphlab::distributed_graph<vertex_data, edge_data> graph_type;


/*
 *The graph loader is used by graph.load to parse lines
 *of the text data file.
 */

bool graph_loader(graph_type& graph, const std::string& fname, const std::string& line ){

	std::stringstream strm(line);

	graphlab::vertex_id_type vid;

	// first entry in the line is a vertex ID
	strm >> vid;

	cout << "line:" << line << endl;
	cout << "vid:" << vid << endl;
	
	float ivalue = 0.0;

	strm >> ivalue;

	graph.add_vertex(vid, vertex_data(vid, ivalue, 0.0, 0.0));

	// while there are elements in the line, continue to read until we fai	
	float vs;
	float vr = 0.0;
	float va = 0.0;
	
	int neighbor_num = 0;
	strm >> neighbor_num;
	//neighbor_num++;
	while(neighbor_num){
		
		neighbor_num--;
		if (strm.fail()) {
			cout << "strm read fail..." << endl;
			break;
		}

		graphlab::vertex_id_type other_vid;
	    	strm >> other_vid;
		strm >> vs;
		
		
		cout << "vid: " << vid 
		     << " other_vid: " << other_vid 
		     << " vs: " << vs 
		     << " vr: " << vr
		     << " va: " << va << endl;
		
		graph.add_edge(vid, other_vid, edge_data(vs, vr, va));

	}
	
	return true;
};


/*
 *This is the message exchange between vertexs
 */
struct message : public graphlab::IS_POD_TYPE {
	
	float s;
	float r;
	float a;

	int source;
	int target;


	message() {}
	explicit message(float fs, float fr, float fa, int sv, int tv): s(fs), r(fr), a(fa), source(sv), target(tv) {}
};


/*
 *This is the gathering type which "accumulates" the a(i, k), s(i, k) and r(,k )
 */

struct set_union_gather : public graphlab::IS_POD_TYPE {

	std::vector<message> messages;
	float maxk;
	/*
	 * Combining with another collection of vertices
	 * union it into the current set,
	 */
	
	set_union_gather& operator+=(const set_union_gather& other) {
		
		for (unsigned int i = 0; i < other.messages.size(); i++) {
			messages.push_back(other.messages[i]);	
			maxk += (other.messages[i].a + other.messages[i].r);
		}
		return *this;
	}
};

class affinity_propagation :
	public graphlab::ivertex_program<graph_type, set_union_gather>,
	public graphlab::IS_POD_TYPE {

private:
	// a variable local to this program
	bool perform_scatter2;
	
public:
	edge_dir_type gather_edges(icontext_type& context, 
				   const vertex_type& vertex) const { 
		/*
		cout << "vertex: " << vertex.data().vertex_id << " OUT_EDGES: " << graphlab::OUT_EDGES << endl;		
		*/
		cout << "graphlab::IN_EDGES" << endl;
		return graphlab::ALL_EDGES;
	} // end of gather_edges 

	gather_type gather(icontext_type& context, const vertex_type& vertex, edge_type& edge) const {

		set_union_gather gather;

		float fs = edge.data().s;
		float fr = edge.data().r;
		float fa = edge.data().a;
				
		int sv = edge.source().data().vertex_id;
		int tv = edge.target().data().vertex_id;

		gather.messages.push_back(message(fs, fr, fa, sv, tv));

		return gather;
	}

	void apply(icontext_type& context, vertex_type& vertex, 
		  const gather_type& total) {
		//cout << "vertex_id: " << vertex.data().vertex_id << endl;

		int center = vertex.data().vertex_id;
		cout << "=====-----=====" << endl;
		cout << "center: " << center << endl;
		
		/*
		 * The target vertexs of current center
		 */
		std::set<int> target_vertexs;
		std::set<int> source_vertexs;

		float maxk = total.maxk;

		// find the target vertexs
		for (unsigned int i = 0; i < total.messages.size(); i++) {
			cout << total.messages[i].s << " "
			     << total.messages[i].r << " "
			     << total.messages[i].a << " "
			     << total.messages[i].source << " " 
			     << total.messages[i].target << endl; 

			if (vertex.data().vertex_id == total.messages[i].source)	
				target_vertexs.insert(total.messages[i].target);

			if (vertex.data().vertex_id == total.messages[i].target)
				source_vertexs.insert(total.messages[i].source);

		}

		former_maxk = maxk;

		for (set<int>::const_iterator iter = target_vertexs.begin(); 
			iter != target_vertexs.end(); 
			++iter) {
		
			float max_as = numeric_limits<float>::min();
			float sum_max_r = 0.0;
			float min_r = numeric_limits<float>::min();

			for (unsigned int i = 0; i < total.messages.size(); i++){
				
				if (total.messages[i].target == center 
					&& total.messages[i].source != *iter){
					
					if (total.messages[i].a + total.messages[i].s < max_as)
						max_as = total.messages[i].a + total.messages[i].s;

					if (total.messages[i].r > 0)
						sum_max_r += total.messages[i].r;
				}					
			
			}

			//TOFIX
			min_r = vertex.data().r + sum_max_r;
	
			vertex.data().msg.push_back(center_msg(*iter, max_as, min_r));			
			vertex.data().a = sum_max_r;
		}

		perform_scatter = false;
		cout << "perform_scatter: " << perform_scatter << endl;
	
	}
		
	edge_dir_type scatter_edges(icontext_type& context,
				    const vertex_type& vertex) const {
		
		if (perform_scatter) {
			cout << "return out edges" << endl;
			return graphlab::OUT_EDGES;
		}
		else {
			cout << "return no edges" << endl;
			return graphlab::NO_EDGES;
		}
		
	}
	
	
	// for each vertex do scatter...
	void scatter(icontext_type& context, const vertex_type& vertex,
		     edge_type& edge) const {
		
		float maxk = 0.0;
		float zero = 0.0;
		for (vector<center_msg>::const_iterator iter = vertex.data().msg.begin(); 
			iter != vertex.data().msg.end(); 
			++iter) {
			if (edge.source().data().vertex_id == vertex.data().vertex_id) {
				if (edge.target().data().vertex_id == (*iter).target_vertex)
					edge.data().r = edge.data().s - (*iter).max_value;
			}
			if (edge.target().data().vertex_id == vertex.data().vertex_id) {
				if (edge.target().data().vertex_id == (*iter).target_vertex)
					edge.data().a = std::min(zero, (*iter).min_value);
			}

			maxk += (edge.data().r + edge.data().a);

		}

		perform_scatter = ((maxk - former_maxk) > 1E-2) ? true : false;
		if (perform_scatter) { 
			cout << "signal edge.target() -------=-------- " << endl;

			context.signal(edge.target()); 	     
		}
	}
	
};

/**
 *We want to save the final graph so we define a write which will be
 *used in graph.save("path/prefix", affinity_propagation_writer()) to save the graph
 */

struct affinity_propagation_writer {

	std::string save_vertex(const graph_type::vertex_type& v) {

		std::stringstream strm;
		//strm << v.data().vertex_id << "\t" << v.data().s << "\n";
		strm << v.data().vertex_id << "\n";

		return strm.str();
	}

	std::string save_edge(graph_type::edge_type e) {
		return "";
	}
};

int main(int argc, char** argv) {
	
	graphlab::mpi_tools::init(argc, argv);
	graphlab::distributed_control dc; 
	
	const std::string description = "Affinity Propagation Clustering";
	graphlab::command_line_options clopts(description);
	std::string exec_type = "sync";
	clopts.attach_option("engine", exec_type, "The type of engine to use {async, sync}.");

	std::string saveprefix;
	clopts.attach_option("saveprefix", saveprefix, "If set, will save the examplar to a sequence of files with prefix saveprefix");

	graph_type graph(dc);
	
	graph.load("ap_graph.txt", graph_loader);

	graph.finalize();

	typedef graphlab::omni_engine<affinity_propagation> engine_type;

	engine_type engine(dc, graph, exec_type, clopts);
	engine.signal_all();
	graphlab::timer timer;
	engine.start();

	// @irwenqiang, because of this line, I spend a total afternoon to
	// figure out the bug! 
	// stupid!..
	// graphlab::mpi_tools::finalize();
	
	const double runtime = timer.current_time();
	dc.cout() 
	<< "-----------------------------------" << endl
	<< "Final Runtime (seconds):	" << runtime
	<< endl
	<< "Updates executed: " << engine.num_updates()
	<< endl
	<< "Update Rate (updates/seconds): "
	<< engine.num_updates() / runtime << endl;

	/*
	if (saveprefix != ""){
		graph.save(saveprefix, affinity_propagation_writer(),
			   false,	// do not gzip
			   true,	// save vertices
			   false);	// do not save edges
	}
	*/
	graphlab::mpi_tools::finalize();

	return EXIT_SUCCESS;
}

