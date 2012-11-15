/*
#include <string>
#include <graphlab.hpp>

struct web_page {
	std::string pagename;
	double pagerank;
	web_page():pagerank(0.0) { }
	explicit web_page(std::string name):pagename(name), pagerank(0.0) { }

	void save(graphlab::oarchive& oarc) const {
		oarc << pagename << pagerank;
	}

	void load(graphlab::iarchive& iarc) {
		iarc >> pagename >> pagerank;
	}
};


typedef graphlab::distributed_graph_type<web_page, graphlab::empty> graph_type;
//typedef graphlab::distributed_graph_type<web_page, graphlab::empty> graph_type;

int main(int argc, char** argv) {
	graphlab::mpi_tools::init(argc, argv);
	graphlab::distributed_control dc;

	dc.cout() << "Hello World!\n";

	graphlab::mpi_tools::finalize();

}
*/

#include <string>
#include <iostream>
#include <graphlab.hpp>

using namespace std;

struct web_page {
	std::string pagename;
	double pagerank;
	web_page():pagerank(0.0) { }
	explicit web_page(std::string name):pagename(name),pagerank(0.0){ }

	void save(graphlab::oarchive& oarc) const {
		oarc << pagename << pagerank;
	}

	void load(graphlab::iarchive& iarc) {
		iarc >> pagename >> pagerank;
	}
};
						 
typedef graphlab::distributed_graph<web_page, graphlab::empty> graph_type;

bool line_parser(graph_type& graph, 
		const std::string& filename,
		const std::string& textline){

	
	std::stringstream strm(textline);

	std::cout << textline << endl;

	graphlab::vertex_id_type vid;

	std::string pagename;

	strm >> vid;
	strm >> pagename;
	cout << "vid:" << vid << endl;
	cout << "pagename:" << pagename << endl;

	graph.add_vertex(vid, web_page(pagename));

	std::string str;

	while(1){
		graphlab::vertex_id_type other_vid;
		strm >> other_vid;
		cout << "other_vid:" << other_vid << endl;
		if (strm.fail()) {
			cout << "strm fail!";
			break;
		}
		graph.add_edge(vid, other_vid);

	}
	return true;
}

int main(int argc, char** argv) {

	graphlab::mpi_tools::init(argc, argv);

	graphlab::distributed_control dc;
							      
	graph_type graph(dc);
	graph.load("graph.txt", line_parser);
	graphlab::mpi_tools::finalize();
}
