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
	#include <boost/program_options.hpp>

	using namespace std;

	// float, optional, default: 0.5
	// Damping factor between 0.5 and 1.
	float damping = 0.5;
	int MAX_ITER = 100;
	string INPUT_FILE = "ap_graph.txt";

	/*
	 * the incoming edge's vertex_id and the max(a(i, k) + s(i, k))
	 */
	struct center_msg : public graphlab::IS_POD_TYPE {
			
		int target_vertex;
		float value;

		center_msg() {}

		explicit center_msg(int a_vertex, float a_value): target_vertex(a_vertex), value(a_value) {}	
	};

	/*
	 * the vertex_data type:
	 * every vertex_data, contains: 
	 * 	1. the max_value:	max(a(i, k) + s(i, k)) 
	 * 	2. the second:	 max(a(i, k) + s(i, k))
	 * 	among all it's incoming edges
	 * 	3. the sum(max(0, r(i, k)))
	 *	among all it's incoming edges
	 *	4. similarity of itself, it is 0.0:	s(k, k)
	 *	5. responsibiliry of itself:		r(k, k)
	 *	6. avaliablity of iteself:		a(k, k)
	 *	7. count:	if count % 2 equals 0, then update the r(i ,k) of vertex_id's outcoming edges
	 *	 		if count % 2 equals 1, then update the a(i, k) of vertex_id's outcoming edges and vertex_id's avaliablity
	 */
	struct vertex_data : public graphlab::IS_POD_TYPE {

		int vertex_id;

		center_msg max_value;
		center_msg sec_max_value;
		
		float sum_r;
		float s;
		float r;
		float a;
		
		int count;

		vertex_data() {}
		explicit vertex_data(int id, float fs, float fr, float fa):vertex_id(id), s(fs), r(fr), a(fa), count(0), sum_r(0.0) { }

	};

	/*
	 * the edge_data type:
	 * similarity between source vertex and target vertex: s
	 * responsibility: r
	 * avaliablity: a
	 */
	struct edge_data : public graphlab::IS_POD_TYPE {

		float s;
		float r;
		float a;

		edge_data(){ }
		explicit edge_data(float vs, float vr, float va):s(vs), r(vr), a(va) { }

	};

	typedef graphlab::distributed_graph<vertex_data, edge_data> graph_type;

	bool graph_loader(graph_type& graph, const std::string& fname, const std::string& line ){

		std::stringstream strm(line);

		graphlab::vertex_id_type vid;

		// first entry in the line is a vertex ID
		strm >> vid;

		float ivalue = 0.0;

		strm >> ivalue;

		graph.add_vertex(vid, vertex_data(vid, ivalue, ivalue, 0.0));

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

			graph.add_edge(other_vid, vid, edge_data(vs, vr, va));

		}
		
		return true;

	};

	/*
	 * the message sum in gather phase
	 * the message store the current vertex's:
	 * 	1. max(a(i, k) + s(i, k)) 
	 * 	2. second max(a(i, k) + s(i, k))
	 * 	3. sum(max(0, r(i, k)))
	 * 	4. count: as that in vertex_data type
	 */
	struct set_union_gather : public graphlab::IS_POD_TYPE {

		int target_vertex;
		float s;
		float r;
		float a;
		
		center_msg max;
		center_msg sec_max;

		float sum;
		int count;
		set_union_gather() {}
		explicit set_union_gather(int a_source, 
					  float a_s, 
					  float a_r, 
					  float a_a,
					  center_msg a_max,
					  center_msg a_sec_max,
					  float a_sum,
					  int a_count): 
			target_vertex(a_source), 
			s(a_s), 
			r(a_r), 
			a(a_a), 
			max(a_max), 
			sec_max(a_sec_max),
			sum(a_sum),
			count(a_count)
			{ }
		
		set_union_gather& operator+=(const set_union_gather& other) {
			
			if (count % 2 == 0){
				float other_as = other.a + other.s;

				if (other_as > max.value) {
					sec_max.value = max.value;
					sec_max.target_vertex = max.target_vertex;

					max.value = other_as;
					max.target_vertex = other.target_vertex;
					
					return *this;
				}
				
				if (other_as <= max.value && other_as >= sec_max.value) {
					sec_max.value = other_as;
					sec_max.target_vertex = other.target_vertex;
				}

				return *this;
			}else {
				sum += std::max(float(0.0), other.r);	
				return *this;
			}
		}
	};

	class affinity_propagation :
		public graphlab::ivertex_program<graph_type, set_union_gather>,
		public graphlab::IS_POD_TYPE {

	public:

		edge_dir_type gather_edges(icontext_type& context, 
					   const vertex_type& vertex) const { 

			// v2: the in edges will be enough
			return graphlab::IN_EDGES;
		}

		/*
		 * compute:
		 * 	1. the max_value:       max(a(i, k) + s(i, k))   
		 * 	2. the sec_max_value:	 max(a(i, k) + s(i, k))    
		 */
		gather_type gather(icontext_type& context, const vertex_type& vertex, edge_type& edge) const {

			float fs = edge.data().s;
			float fr = edge.data().r;
			float fa = edge.data().a;
							
			int sv = edge.source().data().vertex_id;
			
			return set_union_gather(sv, fs, fr, fa, center_msg(sv, edge.data().s + edge.data().a), center_msg(sv, -numeric_limits<float>::max()), std::max(float(0.0), fr),vertex.data().count);
		
		}

		void apply(icontext_type& context, vertex_type& vertex, 
			  const gather_type& total) {
			
			if (vertex.data().count % 2 == 0) {
				//cout << "apply 2" << endl;
				vertex.data().max_value.target_vertex = total.max.target_vertex;
				vertex.data().max_value.value = total.max.value;

				vertex.data().sec_max_value.target_vertex = total.sec_max.target_vertex;
				vertex.data().sec_max_value.value = total.sec_max.value;
			}else {
				//cout << "apply 1" << endl;
				vertex.data().sum_r = total.sum;
				vertex.data().a = total.sum;
			}

			vertex.data().count += 1;
		
		}
			
		edge_dir_type scatter_edges(icontext_type& context,
					    const vertex_type& vertex) const {
			
			return graphlab::OUT_EDGES;		
		}
		
		void scatter(icontext_type& context, const vertex_type& vertex,
			     edge_type& edge) const {
					
			if ((vertex.data().count - 1) % 2 == 0){
				float old_r = edge.data().r;
				if (edge.target().data().vertex_id != vertex.data().max_value.target_vertex) {
					edge.data().r = vertex.data().r - vertex.data().max_value.value;
				}
				else{
					edge.data().r = vertex.data().r - vertex.data().sec_max_value.value;
				}
			
				edge_data().r = (1 - damping) * edge_data().r + damping * old_r;
			}else {
				float old_a = edge.data().a;

				edge.data().a = std::min(float(0.0), vertex.data().r + vertex.data().sum_r - std::max(float(0.0),edge.data().r));
				edge.data().a = (1 - damping) * edge.data().a + damping * old_a;
			}
			
			if (vertex.data().count < MAX_ITER) {

				context.signal(edge.target());

				if (vertex.data().a + vertex.data().r > 0)
					cout << "examplar: " << vertex.data().vertex_id << endl;
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
	std::string exec_type = "async";
	clopts.attach_option("engine", exec_type, "The type of engine to use {async, sync}.");

	std::string saveprefix;
	clopts.attach_option("saveprefix", saveprefix, "If set, will save the examplar to a sequence of files with prefix saveprefix");

	std::string tolerance;
	clopts.attach_option("tolerance", tolerance, "the tolerance to terminate the computation");
	
	clopts.attach_option("inputfile", INPUT_FILE, "the input file of the graph");
	clopts.attach_option("max_iter", MAX_ITER, "max iterations");
	
	if(!clopts.parse(argc, argv)) {

		dc.cout() << "Error in parsing command line arguments." << std::endl;
		return EXIT_FAILURE;
	}

	graph_type graph(dc, clopts);

	graph.load(INPUT_FILE, graph_loader);

	graph.finalize();

	typedef graphlab::omni_engine<affinity_propagation> engine_type;
	
	engine_type engine(dc, graph, exec_type, clopts);
	engine.signal_all();
	graphlab::timer timer;
	engine.start();
	
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

