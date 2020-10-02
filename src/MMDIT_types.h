#ifndef GRAPHLINE_TYPES_H
#define GRAPHLINE_TYPES_H 1

#include <set>
#include <string>
#include <vector>
#include <stdint.h>

/**A node in a sequence graph (with epsilon jumps).*/
typedef struct{
	/**The number of characters at this node: zero for an epsilon node.*/
	uintptr_t numChars;
	/**The caracters for this node. NOT NULL TERMINATED.*/
	const char* nodeChars;
	/**The number of back links this has.*/
	uintptr_t numBack;
	/**The back link indices.*/
	uintptr_t* backLinks;
	/**The number of forward links this has.*/
	uintptr_t numForw;
	/**The forward link indices.*/
	uintptr_t* forwLinks;
} GraphlineEJNode;

/**A sequence graph.*/
class GraphlineEJGraph{
public:

	/**Set up an empty graph.*/
	GraphlineEJGraph();
	/**Tear down a graph.*/
	~GraphlineEJGraph();

	/**
	 * Add a node to the graph.
	 * @param begEndFlags The flags for whether the node can begin or end.
	 * @param nodeSeq The sequence for the node.
	 * @param nodeBacks The back links for the node.
	 * @param refS The index in the reference this corresponds to: -1 for stuff not in the reference.
	 * @param refE The index in the reference after this sequence: -1 for stuff not in the reference.
	 */
	void addNode(int begEndFlags, std::string* nodeSeq, std::vector<uintptr_t>* nodeBacks, uintptr_t refS, uintptr_t refE);

	/**
	 * Take all the disjointed node data and make it useful.
	 */
	void compileNodes();

	/**The nodes in this graph.*/
	std::vector<GraphlineEJNode> allNodes;
	/**The nodes that can be at the end.*/
	std::set<uintptr_t> endNodes;
	/**The nodes at which this can begin.*/
	std::set<uintptr_t> beginNodes;

	/**Storage for sequence data for the nodes.*/
	std::string storeSeqs;
	/**Storage for the forward data for the nodes.*/
	std::vector<uintptr_t> storeForws;
	/**Storage for the backward data for the nodes.*/
	std::vector<uintptr_t> storeBacks;

	/**The location in the linear reference each node corresponds to.*/
	std::vector< std::pair<uintptr_t,uintptr_t> > nodeRefCor;
};

/**An entry in a graph traversal.*/
typedef struct{
	/**The next index in the sequence to handle.*/
	uintptr_t seqInd;
	/**The node to handle next.*/
	uintptr_t nodeInd;
	/**The index in that node to handle next.*/
	uintptr_t nodeSeqInd;
	/**The previous entry for this traversal, one indexed. Zero for none.*/
	uintptr_t prevEnt;
	/**The number of errors to get to this point.*/
	uintptr_t offCost;
} GraphlineTraverseEntry;

/**A storage for a traversal.*/
class GraphlineTraversalDump{
public:
	/**The graph used in building this thing.*/
	uintptr_t graphAddr;
	/**The edit distance used during the traversal.*/
	int maxEdit;
	/**The traversal options.*/
	std::vector<GraphlineTraverseEntry> travGraph;
	/**The options that represent full traversals.*/
	std::vector<uintptr_t> travEnds;
	/**The original query.*/
	std::string origQuery;
};

#endif
