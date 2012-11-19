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
float tolerance = 10;
// float, optional, default: 0.5
// Damping factor between 0.5 and 1.
float damping = 0.5;
bool isMax = false;

struct center_msg : public graphlab::IS_POD_TYPE {

	int target_vertex;
	//max(a(i, k') + s(i, k')),	s.t. k'!= k
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
	int count;

	vertex_data() { }
	explicit vertex_data(int id, float fs, float fr, float fa):vertex_id(id), s(fs), r(fr), a(fa), count(0) { }

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

	float ivalue = 0.0;

	strm >> ivalue;

	graph.add_vertex(vid, vertex_data(vid, ivalue, 0.0, 0.0));

	float vs;
	float vr = 0.0;
	float va = 0.0;
	
	int neighbor_num = 0;

	strm >> neighbor_num;

	// while there are elements in the line, continue to read until we fail	
	while(neighbor_num){
		
		neighbor_num--;
		if (strm.fail()) {
			cout << "strm read fail..." << endl;
			break;
		}

		graphlab::vertex_id_type other_vid;
	    	strm >> other_vid;
		strm >> vs;

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

	// accumulate the what? @irwenqiang
	float maxk;
	/*
	 * Override the +=/sum operator in the pharse between 
	 * gather and apply:
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

public:
	edge_dir_type gather_edges(icontext_type& context, 
				   const vertex_type& vertex) const { 

		// currently, just consider handling the all edges,
		// the optimization will be taked in to account lately.
		// @irwenqiang

		// v2: the in edges will be enough
		return graphlab::IN_EDGES;
	}

	/**
	 * gather all the edges that have connnections to the current vertex,i.e. 
	 * the center
	 */
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

	/**
	 * According the signal of the apply function,
	 * this apply pharse operates on the the current vertex,
	 * that is operating on the center
	 */
	void apply(icontext_type& context, vertex_type& vertex, 
		  const gather_type& total) {

		vertex.data().count++;

		if (vertex.data().count > 10)
			isMax = true;

		int center = vertex.data().vertex_id;
				
		/*
		 * The target vertexs of current center
		 */
		std::set<int> in_vertexs;

		float maxk = total.maxk;

		// find the target vertexs
		for (unsigned int i = 0; i < total.messages.size(); i++) {

			if (center == total.messages[i].target)
				in_vertexs.insert(total.messages[i].source);
		}

		former_maxk = maxk;

		/**
		 * the ground truth is that 
		 * in the affinity propagation, firstly, consider all the link input
		 * edges, and then, update the value of the link output edges
		 */
		for (set<int>::const_iterator iter = in_vertexs.begin(); 
			iter != in_vertexs.end(); 
			++iter) {
				
			float max_as = -(numeric_limits<float>::max() - 2);
			float sum_max_r = 0.0;
			float min_r = -(numeric_limits<float>::max() - 2);
			float zero = 0.0;
			for (unsigned int i = 0; i < total.messages.size(); i++){
				
				if (total.messages[i].source != *iter){
					
					if (total.messages[i].a + total.messages[i].s > max_as)
						max_as = total.messages[i].a + total.messages[i].s;

					if (total.messages[i].r > 0)
						sum_max_r += total.messages[i].r;
				}					
			
			}
			
			cout << "sum_max_r: " << sum_max_r << endl;
			min_r  = std::min(vertex.data().r + sum_max_r, zero);
			

			// update the value of the vertex data in this apply pharse
			// *iter will be the k of r(i, k)
			vertex.data().msg.push_back(center_msg(*iter, max_as, min_r));			
			// update the avaliablity of the vertex, according to
			// a(k, k) = sum(max(0, r(i', k))), subject to i' != k
			vertex.data().a = sum_max_r;

			// there's not any infomation about how to update the responsibility of the vertex k
			// so, vertex.data().r = ?? @irwenqiang
		}
	}
		
	// what's the scatter_edges done?
	// or
	// what's the default operation of the scatter_edges pharse?
	// @irwenqiang
	
	edge_dir_type scatter_edges(icontext_type& context,
				    const vertex_type& vertex) const {
		
		return graphlab::OUT_EDGES;		
	}
	
	// It's time to update the edge_data in the scatter pharse...
	void scatter(icontext_type& context, const vertex_type& vertex,
		     edge_type& edge) const {
		
		float maxk = 0.0;
		float zero = 0.0;

		float old_r = edge.data().r;
		float old_a = edge.data().a;

		// in the scatter pharse, due to arrange in the gather and, specially 
		// in the apply pharse, we can easily update the edge datas of the center's
		// link output edges
		for (vector<center_msg>::const_iterator outer_iter = vertex.data().msg.begin(); 
			outer_iter != vertex.data().msg.end(); 
			++outer_iter) {
			
			if (edge.target().data().vertex_id == (*outer_iter).target_vertex) {
				edge.data().r = edge.data().s - (*outer_iter).max_value;
				cout << "edge.data().s: " << edge.data().s << "\t"
				     << "max: " << (*outer_iter).max_value
				     << endl;
				edge.data().r = (1 - damping) * edge.data().r + damping * old_r;
				edge.data().a = std::min(zero, (*outer_iter).min_value);
				edge.data().a = (1 - damping) * edge.data().a + damping * old_a;
			}

			maxk += (edge.data().r + edge.data().a);

		}
//		perform_scatter = (std::fabs(maxk - former_maxk) < 1E-6);
		perform_scatter = (std::fabs(std::fabs(maxk) - std::fabs(former_maxk)) < tolerance);
		
		if (isMax)
			return;
		if (!perform_scatter) { 
			context.signal(edge.target()); 	     
		}
		
		if (perform_scatter) {
			
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
		strm << v.data().vertex_id 
		     << "\t"
		     << v.data().s
		     << "\t"
		     << v.data().r
		     << "\t"
		     << v.data().a
		     << endl;

		return strm.str();
	}

	std::string save_edge(graph_type::edge_type e) {
		std::stringstream strm;
		strm << e.source().data().vertex_id
		     << "\t"
		     << e.target().data().vertex_id
		     << "\t"
		     << e.data().s 
		     << "\t"
		     << e.data().r
		     << "\t"
		     << e.data().a
		     << endl;

		return strm.str();

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

	std::string tolerance;
	clopts.attach_option("tolerance", tolerance, "the tolerance to terminate the computation");
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

	
	saveprefix = "output";
	if (saveprefix != ""){
		cout << "savvvvvvvvvvvvvvvvvvvvvvvvvvveprefix" << endl;
		graph.save(saveprefix, affinity_propagation_writer(),
			   false,	// do not gzip
			   false,	// save vertices
			   true);	// do not save edges
	}

	graphlab::mpi_tools::finalize();

	return EXIT_SUCCESS;
}

