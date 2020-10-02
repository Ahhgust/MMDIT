
#include <iostream>
#include <algorithm>

#include "MMDIT_types.h"

#include <set>
#include <string>
#include <vector>
#include <stdint.h>



#define GRAPHLINE_NODE_BEGIN 1
#define GRAPHLINE_NODE_END 2

GraphlineEJGraph::GraphlineEJGraph(){
	//nothing to do
}

GraphlineEJGraph::~GraphlineEJGraph(){
	//nothing to do
}

void GraphlineEJGraph::addNode(int begEndFlags, std::string* nodeSeq, std::vector<uintptr_t>* nodeBacks, uintptr_t refS, uintptr_t refE){
	uintptr_t thisNodeInd = allNodes.size();
	//pointers are set to integers (the offsets in the array), will correct in compileNodes.
	//also do not have forward links, will also correct later
	GraphlineEJNode newNode = {nodeSeq->size(), (const char*)storeSeqs.size(), nodeBacks->size(), (uintptr_t*)storeBacks.size(), 0, 0};
	allNodes.push_back(newNode);
	if(begEndFlags & GRAPHLINE_NODE_BEGIN){
		beginNodes.insert(thisNodeInd);
	}
	if(begEndFlags & GRAPHLINE_NODE_END){
		endNodes.insert(thisNodeInd);
	}
	storeSeqs.append(*nodeSeq);
	storeBacks.insert(storeBacks.end(), nodeBacks->begin(), nodeBacks->end());
	nodeRefCor.push_back( std::pair<uintptr_t,uintptr_t>(refS, refE) );
}

void GraphlineEJGraph::compileNodes(){
	//fix the pointers
	for(uintptr_t i = 0; i<allNodes.size(); i++){
		GraphlineEJNode* curNode = &(allNodes[i]);
		curNode->nodeChars = &(storeSeqs[(uintptr_t)(curNode->nodeChars)]);
		curNode->backLinks = &(storeBacks[(uintptr_t)(curNode->backLinks)]);
	}
	//build the forward links
	std::vector< std::pair<uintptr_t,uintptr_t> > linkStore;
	for(uintptr_t i = 0; i<allNodes.size(); i++){
		GraphlineEJNode* curNode = &(allNodes[i]);
		for(uintptr_t j = 0; j<curNode->numBack; j++){
			linkStore.push_back( std::pair<uintptr_t,uintptr_t>(curNode->backLinks[j], i) );
		}
	}
	std::sort(linkStore.begin(), linkStore.end());
	storeForws.resize(linkStore.size());
	for(uintptr_t i = 0; i<linkStore.size(); i++){
		storeForws[i] = linkStore[i].second;
	}
	uintptr_t curForwOff = 0;
	for(uintptr_t i = 0; i<allNodes.size(); i++){
		GraphlineEJNode* curNode = &(allNodes[i]);
		curNode->forwLinks = &(storeForws[curForwOff]);
		uintptr_t endForwOff = curForwOff;
		while((endForwOff < linkStore.size()) && (linkStore[endForwOff].first == i)){
			endForwOff++;
		}
		curNode->numForw = (endForwOff - curForwOff);
		curForwOff = endForwOff;
	}

	/*
	for(uintptr_t i = 0; i<allNodes.size(); i++){
		GraphlineEJNode* curNode = &(allNodes[i]);
		std::string toOut(curNode->nodeChars, curNode->nodeChars + curNode->numChars);
		std::cout << i << "\t" << '"' << toOut << '"';
		for(uintptr_t j = 0; j<curNode->numForw; j++){
			std::cout << " " << curNode->forwLinks[j];
		}
		std::cout << std::endl;
	}
	std::cout << "B ";
	for(std::set<uintptr_t>::iterator setIt = beginNodes.begin(); setIt != beginNodes.end(); setIt++){
		std::cout << " " << *setIt;
	}
	std::cout << std::endl;
	std::cout << "E ";
	for(std::set<uintptr_t>::iterator setIt = endNodes.begin(); setIt != endNodes.end(); setIt++){
		std::cout << " " << *setIt;
	}
	std::cout << std::endl;
	*/
}

/**
 * This will traverse a sequence graph with a sequence, finding the matching traversals within some distance of perfect.
 * @param maxEdit The maximum edit distance.
 * @param forSeq The sequence to walk down the graph with.
 * @param forGraph The graph to walk down.
 * @param travGraph The traversals through the graph.
 * @param possibleEnds The traversal entry indices that meet the criteria.
 */
void graphlineTraverseSequenceGraph(uintptr_t maxEdit, std::string* forSeq, GraphlineEJGraph* forGraph, std::vector<GraphlineTraverseEntry>* travGraph, std::set<uintptr_t>* possibleEnds){
	std::vector<uintptr_t> curHor;
	std::vector<uintptr_t> nextHor;
	std::vector<uintptr_t> tmpHor;
	//put all the starting pieces on the horizon
	for(std::set<uintptr_t>::iterator begIt = forGraph->beginNodes.begin(); begIt!=forGraph->beginNodes.end(); begIt++){
		GraphlineTraverseEntry curEnt = {0,*begIt,0,0,0};
		curHor.push_back(travGraph->size());
		travGraph->push_back(curEnt);
	}
	//run down the sequence
	for(uintptr_t si = 0; si<=forSeq->size(); si++){
		//expand epsilon jumps and end of nodes
		tmpHor.clear();
		tmpHor.insert(tmpHor.end(), curHor.begin(), curHor.end());
		uintptr_t thI = 0;
		while(thI < tmpHor.size()){
			GraphlineTraverseEntry horEnt = ((*travGraph)[tmpHor[thI]]);
			GraphlineEJNode* horNode = &(forGraph->allNodes[horEnt.nodeInd]);
			if(horEnt.nodeSeqInd >= horNode->numChars){
				for(uintptr_t i = 0; i<horNode->numForw; i++){
					uintptr_t neNode = horNode->forwLinks[i];
					uintptr_t nePrev = tmpHor[thI];
					uintptr_t neCost = horEnt.offCost;
					GraphlineTraverseEntry nextEnt = {si, neNode, 0, nePrev+1, neCost};
					tmpHor.push_back(travGraph->size());
					travGraph->push_back(nextEnt);
				}
			}
			else{
				//consider skipping the reference
				if(horEnt.offCost + 1 <= maxEdit){
					GraphlineTraverseEntry nextEnt = {si, horEnt.nodeInd, horEnt.nodeSeqInd+1, tmpHor[thI]+1, horEnt.offCost+1};
					tmpHor.push_back(travGraph->size());
					travGraph->push_back(nextEnt);
				}
			}
			thI++;
		}
		if(si < forSeq->size()){
			//try to follow the sequence character
			char curChar = (*forSeq)[si];
			nextHor.clear();
			for(thI = 0; thI<tmpHor.size(); thI++){
				GraphlineTraverseEntry horEnt = ((*travGraph)[tmpHor[thI]]);
				GraphlineEJNode* horNode = &(forGraph->allNodes[horEnt.nodeInd]);
				if(horEnt.nodeSeqInd >= horNode->numChars){ continue; }
				//look for match/mismatch
				GraphlineTraverseEntry nextEnt = {si+1, horEnt.nodeInd, horEnt.nodeSeqInd+1, tmpHor[thI]+1, horEnt.offCost};
				if(curChar != horNode->nodeChars[horEnt.nodeSeqInd]){
					nextEnt.offCost++;
				}
				if(nextEnt.offCost <= maxEdit){
					nextHor.push_back(travGraph->size());
					travGraph->push_back(nextEnt);
				}
				//consider skipping the sequence
				if(horEnt.offCost + 1 <= maxEdit){
					GraphlineTraverseEntry skipEnt = {si+1, horEnt.nodeInd, horEnt.nodeSeqInd, tmpHor[thI]+1, horEnt.offCost+1};
					nextHor.push_back(travGraph->size());
					travGraph->push_back(skipEnt);
				}
			}
		}
		else{
			std::swap(tmpHor,nextHor);
		}
		//prepare for the next item
		std::swap(curHor,nextHor);
	}
	//look for everything that ends at the end of an ending node
	for(uintptr_t i = 0; i<curHor.size(); i++){
		//std::cout << curHor[i] << std::endl;
		GraphlineTraverseEntry* horEnt = &((*travGraph)[curHor[i]]);
		GraphlineEJNode* horNode = &(forGraph->allNodes[horEnt->nodeInd]);
		if(forGraph->endNodes.count(horEnt->nodeInd) == 0){continue;}
		if(horEnt->nodeSeqInd < horNode->numChars){continue;}
		possibleEnds->insert(curHor[i]);
	}

	/*
	std::cout << "Traversal " << *forSeq << std::endl;
	for(uintptr_t i = 0; i<travGraph->size(); i++){
		GraphlineTraverseEntry curEnt = (*travGraph)[i];
		std::cout << i << " " << curEnt.seqInd << " " << curEnt.nodeInd << " " << curEnt.nodeSeqInd << " " << (curEnt.prevEnt-1) << " " << curEnt.offCost << std::endl;
	}
	for(std::set<uintptr_t>::iterator setIt = possibleEnds->begin(); setIt != possibleEnds->end(); setIt++){
		std::cout << *setIt << std::endl;
	}
	*/
}

/**
 * Note which nodes and edges are missing from the graph.
 * @param forGraph The graph being traversed.
 * @param travGraph The traversal nodes.
 * @param travInd The index of the terminal traversal node.
 * @param hitNodes The nodes touched by the traversal.
 * @param hitEdges The edges touched by the traversal.
 */
void graphlineNoteMissedNodeEdge(GraphlineEJGraph* forGraph, std::vector<GraphlineTraverseEntry>* travGraph, uintptr_t travInd, std::set<uintptr_t>* hitNodes, std::set<uintptr_t>* hitEdges){
	uintptr_t curEntI = travInd + 1;
	while(curEntI){
		curEntI = curEntI - 1;
		GraphlineTraverseEntry* curEnt = &((*travGraph)[curEntI]);
		uintptr_t curNodeI = curEnt->nodeInd;
		hitNodes->insert(curNodeI);
		if(curEnt->prevEnt){
			uintptr_t prevNodeI = (*travGraph)[curEnt->prevEnt].nodeInd;
			if(prevNodeI != curNodeI){
				GraphlineEJNode* prevNode = &(forGraph->allNodes[prevNodeI]);
				uintptr_t* selLink = std::lower_bound(prevNode->forwLinks, prevNode->forwLinks+prevNode->numForw, curNodeI);
				hitEdges->insert(selLink - &(forGraph->storeForws[0]));
			}
		}
		curEntI = curEnt->prevEnt;
	}
}

//*****************************************************************************
//RCpp interface

#include <Rcpp.h>

//' Given a description of the variants in a mixture, will generate a graph for use with later analysis.
//'
//' @param refRS The reference sequence.
//' @param varS The starting points of the variant regions.
//' @param varE The ending points of the variant regions.
//' @param varSeq The sequences for each variant.
//' @return An opaque object containing the graph.
// [[Rcpp::export]]
Rcpp::XPtr<GraphlineEJGraph> makeVariantGraph(Rcpp::String refRS, Rcpp::IntegerVector varS, Rcpp::IntegerVector varE, Rcpp::StringVector varSeq){
	if(varS.size() != varE.size()){
		throw std::domain_error("Variant vectors must have the same length.");
	}
	if(varS.size() != varSeq.size()){
		throw std::domain_error("Variant vectors must have the same length.");
	}
	uintptr_t refLen = strlen(refRS.get_cstring());
	//set up the ranges
	std::map< std::pair<uintptr_t,uintptr_t>, std::vector<std::string> > regVars;
	std::vector< std::pair<uintptr_t,uintptr_t> > allRegs;
	uintptr_t numV = varS.size();
	for(uintptr_t i = 0; i<numV; i++){
		if(varS[i] < 0){ throw std::domain_error("Cannot have negative indices in the reference."); }
		if(varE[i] < 0){ throw std::domain_error("Cannot have negative indices in the reference."); }
		if(varS[i] > refLen){ throw std::domain_error("Cannot have variants beyond the end of the reference."); }
		if(varE[i] > refLen){ throw std::domain_error("Cannot have variants beyond the end of the reference."); }
		if(varE[i] < varS[i]){ throw std::domain_error("Variants cannot end before they begin."); }
		std::pair<uintptr_t,uintptr_t> varRng(varS[i], varE[i]);
		std::string seqAStr(varSeq[i]);
		regVars[varRng].push_back(seqAStr);
		allRegs.push_back(varRng);
	}
	std::sort(allRegs.begin(), allRegs.end());
	//make sure the regions are disjoint
	std::set< std::pair<uintptr_t,uintptr_t> > uniqRegs;
	for(uintptr_t i = 0; i<numV; i++){
		std::pair<uintptr_t,uintptr_t> regCur = allRegs[i];
		if(i){
			std::pair<uintptr_t,uintptr_t> regPre = allRegs[i-1];
			if(regCur.first < regPre.second){
				if((regCur.first != regPre.first) || (regCur.second != regPre.second)){
					throw std::domain_error("Variants with overlapping regions.");
				}
			}
		}
		uniqRegs.insert(regCur);
	}
	//build up the graph
	GraphlineEJGraph* endGraph = new GraphlineEJGraph();
	std::string tmpSeq;
	std::vector<uintptr_t> tmpBack;
	endGraph->addNode(GRAPHLINE_NODE_BEGIN, &tmpSeq, &tmpBack, 0, 0);
	uintptr_t lastNodeI = 0;
	uintptr_t curRefP = 0;
	for(std::set< std::pair<uintptr_t,uintptr_t> >::iterator setIt = uniqRegs.begin(); setIt != uniqRegs.end(); setIt++){
		std::pair<uintptr_t,uintptr_t> curReg = *setIt;
		if((curRefP < curReg.first) && (curRefP < refLen)){
			tmpSeq.append(refRS, curRefP, curReg.first - curRefP);
			tmpBack.push_back(lastNodeI);
			endGraph->addNode(0, &tmpSeq, &tmpBack, curRefP, curReg.first);
			tmpSeq.clear(); tmpBack.clear(); lastNodeI++;
		}
		tmpBack.push_back(lastNodeI);
		std::vector< std::string >* curVars = &(regVars[curReg]);
		for(uintptr_t i = 0; i<curVars->size(); i++){
			endGraph->addNode(0, &((*curVars)[i]), &tmpBack, curReg.first, curReg.second);
		}
		tmpBack.clear();
		for(uintptr_t i = 0; i<curVars->size(); i++){
			tmpBack.push_back(lastNodeI + i + 1);
		}
		endGraph->addNode(0, &tmpSeq, &tmpBack, curReg.second, curReg.second);
		tmpBack.clear();
		lastNodeI += curVars->size() + 1;
		curRefP = curReg.second;
	}
	if(curRefP < refLen){
		tmpSeq.append(refRS, curRefP, refLen - curRefP);
		tmpBack.push_back(lastNodeI);
		endGraph->addNode(0, &tmpSeq, &tmpBack, curRefP, refLen);
		tmpSeq.clear(); tmpBack.clear(); lastNodeI++;
	}
	tmpBack.push_back(lastNodeI);
	endGraph->addNode(GRAPHLINE_NODE_END, &tmpSeq, &tmpBack, refLen, refLen);
	//return
	endGraph->compileNodes();
	Rcpp::XPtr<GraphlineEJGraph> toRet(endGraph);
	return toRet;
}

//' Will return all alignments of a sequence to a graph (within some maximum edit distance).
//'
//' @param sgrap The graph.
//' @param query The sequence.
//' @param maxEdit The maximum edit distance.
//' @return An opaque object containing the graph traversals.
// [[Rcpp::export]]
Rcpp::XPtr<GraphlineTraversalDump> traverseSequenceGraph(Rcpp::XPtr<GraphlineEJGraph> sgrap, Rcpp::String query, int maxEdit){
	GraphlineTraversalDump* travS = new GraphlineTraversalDump();
	travS->maxEdit = maxEdit;
	travS->graphAddr = (uintptr_t)(sgrap.get());
	travS->origQuery = std::string(query.get_cstring());
	std::set<uintptr_t> travEnds;
	std::string queryAS = query;
	graphlineTraverseSequenceGraph(maxEdit, &queryAS, sgrap.get(), &(travS->travGraph), &travEnds);
	travS->travEnds.insert(travS->travEnds.end(), travEnds.begin(), travEnds.end());
	Rcpp::XPtr<GraphlineTraversalDump> toRet(travS);
	return toRet;
}

//' Vectorized version of traverseSequenceGraph.
//'
//' @param sgrap The graph.
//' @param queries The sequences.
//' @param maxEdit The maximum edit distance.
//' @return A list of opaque objects containing the graph traversals.
// [[Rcpp::export]]
Rcpp::List traverseSequencesGraph(Rcpp::XPtr<GraphlineEJGraph> sgrap, Rcpp::StringVector queries, int maxEdit){
	Rcpp::List toRet(queries.size());
	for(uintptr_t i = 0; i<queries.size(); i++){
		Rcpp::String query = queries[i];
		toRet[i] = traverseSequenceGraph(sgrap, query, maxEdit);
	}
	return toRet;
}

//' Will return the edit distances for the possible alignments to a graph.
//'
//' @param toEx The graph traversals.
//' @return The edit distance of each traversal: length zero means no valid traversals were possible.
// [[Rcpp::export]]
Rcpp::IntegerVector getTraversalEditDistances(Rcpp::XPtr<GraphlineTraversalDump> toEx){
	Rcpp::IntegerVector toRet;
	GraphlineTraversalDump* toExP = toEx.get();
	for(uintptr_t i = 0; i<toExP->travEnds.size(); i++){
		toRet.push_back(toExP->travGraph[toExP->travEnds[i]].offCost);
	}
	return toRet;
}

//' Vectorized version of getTraversalEditDistances.
//'
//' @param toExS The graph traversals.
//' @return The edit distance of each traversal: length zero means no valid traversals were possible.
// [[Rcpp::export]]
Rcpp::List getTraversalsEditDistances(Rcpp::List toExS){
	Rcpp::List toRet(toExS.size());
	for(uintptr_t i = 0; i<toExS.size(); i++){
		Rcpp::XPtr<GraphlineTraversalDump> toEx = Rcpp::as< Rcpp::XPtr<GraphlineTraversalDump> >(toExS[i]);
		toRet[i] = getTraversalEditDistances(toEx);
	}
	return toRet;
}

//' Will return the nodes in the graph that the selected traversal hit..
//'
//' @param sgrap The graph.
//' @param toEx The graph traversals.
//' @param travInd The index of the traversal to examine (one-indexed).
//' @return The graph nodes hit by the traversal.
// [[Rcpp::export]]
Rcpp::IntegerVector getTraversalHitNodes(Rcpp::XPtr<GraphlineEJGraph> sgrap, Rcpp::XPtr<GraphlineTraversalDump> toEx, int travInd){
	GraphlineTraversalDump* toExP = toEx.get();
	GraphlineEJGraph* sgrapP = sgrap.get();
	if(toExP->graphAddr != (uintptr_t)sgrapP){ throw std::domain_error("Traversal did not come from the given graph."); }
	if(travInd <= 0){ throw std::domain_error("Attempt to get index <= 0."); }
	if(travInd > toExP->travEnds.size()){ throw std::domain_error("Attempt to get traversal that does not exist."); }
	std::set<uintptr_t> hitNodes;
	std::set<uintptr_t> hitEdges;
	graphlineNoteMissedNodeEdge(sgrapP, &(toExP->travGraph), toExP->travEnds[travInd-1], &hitNodes, &hitEdges);
	Rcpp::IntegerVector toRet;
	for(std::set<uintptr_t>::iterator setIt = hitNodes.begin(); setIt != hitNodes.end(); setIt++){
		toRet.push_back(*setIt);
	}
	return toRet;
}

//' Will return the nodes in the graph that the selected traversals collectively hit.
//'
//' @param sgrap The graph.
//' @param toExS The graph traversals.
//' @param travInds The indices for each traversal.
//' @return The graph nodes hit by the traversals.
// [[Rcpp::export]]
Rcpp::IntegerVector getTraversalsHitNodes(Rcpp::XPtr<GraphlineEJGraph> sgrap, Rcpp::List toExS, Rcpp::IntegerVector travInds){
	if(toExS.size() != travInds.size()){ throw std::domain_error("Traversal indices do not match the traversal list."); }
	GraphlineEJGraph* sgrapP = sgrap.get();
	std::set<uintptr_t> hitNodes;
	std::set<uintptr_t> hitEdges;
	for(uintptr_t i = 0; i<toExS.size(); i++){
		int travInd = travInds[i];
		Rcpp::XPtr<GraphlineTraversalDump> toEx = Rcpp::as< Rcpp::XPtr<GraphlineTraversalDump> >(toExS[i]);
		GraphlineTraversalDump* toExP = toEx.get();
		if(toExP->graphAddr != (uintptr_t)sgrapP){ throw std::domain_error("Traversal did not come from the given graph."); }
		if(travInd <= 0){ throw std::domain_error("Attempt to get index <= 0."); }
		if(travInd > toExP->travEnds.size()){ throw std::domain_error("Attempt to get traversal that does not exist."); }
		graphlineNoteMissedNodeEdge(sgrapP, &(toExP->travGraph), toExP->travEnds[travInd-1], &hitNodes, &hitEdges);
	}
	Rcpp::IntegerVector toRet;
	for(std::set<uintptr_t>::iterator setIt = hitNodes.begin(); setIt != hitNodes.end(); setIt++){
		toRet.push_back(*setIt);
	}
	return toRet;
}

//' Lists the nodes in a sequence graph.
//'
//' @param sgrap The graph.
//' @return The graph node indices.
// [[Rcpp::export]]
Rcpp::IntegerVector listGraphNodes(Rcpp::XPtr<GraphlineEJGraph> sgrap){
	GraphlineEJGraph* sgrapP = sgrap.get();
	Rcpp::IntegerVector toRet;
	for(uintptr_t i = 0; i<sgrapP->allNodes.size(); i++){
		toRet.push_back(i);
	}
	return toRet;
}

//' Enumerate the possible sets of traversals.
//'
//' @param toExS The graph traversals.
//' @return The possible collective traversal indices.
// [[Rcpp::export]]
Rcpp::List enumerateTraversalPossibilities(Rcpp::List toExS){
	//count how many there are
	std::vector<int> maxInds;
	uintptr_t totPos = 1;
	for(uintptr_t i = 0; i<toExS.size(); i++){
		Rcpp::XPtr<GraphlineTraversalDump> toEx = Rcpp::as< Rcpp::XPtr<GraphlineTraversalDump> >(toExS[i]);
		GraphlineTraversalDump* toExP = toEx.get();
		maxInds.push_back(toExP->travEnds.size());
		totPos = totPos * toExP->travEnds.size();
	}
	//build the list
	Rcpp::List toRet(totPos);
	std::vector<int> curInds;
		curInds.insert(curInds.end(), toExS.size(), 1);
	for(uintptr_t i = 0; i<totPos; i++){
		Rcpp::IntegerVector curPos = Rcpp::wrap(curInds);
		toRet[i] = curPos;
		int fi = 0;
		while(fi < curInds.size()){
			curInds[fi]++;
			if(curInds[fi] <= maxInds[fi]){ break; }
			curInds[fi] = 0;
			fi++;
		}
	}
	return toRet;
}

/**
 * Helper method: avoids redundant calculation.
 * @param sgrapP THe graph.
 * @param toExS The graph traversals.
 * @param travISet The sets of traversals.
 * @param allCount The place to put the results.
 * @param fromQ The query index.
 * @param minTI The smallest traversal set index to consider.
 * @param maxTI The index to stop at.
 * @param hitSoFar The nodes that have been hit up to the fromQ query.
 */
void countTraversalsMissedNodesHelp(GraphlineEJGraph* sgrapP, Rcpp::List toExS, Rcpp::List travISet, std::vector<int>* allCount, int fromQ, int minTI, int maxTI, std::set<uintptr_t>* hitSoFar){
	if(fromQ >= toExS.size()){
		int nodeMiss = sgrapP->allNodes.size() - hitSoFar->size();
		for(int i = minTI; i<maxTI; i++){
			allCount->push_back(nodeMiss);
		}
		return;
	}
	int curI = minTI;
	while(curI < maxTI){
		std::set<uintptr_t> curHits = *hitSoFar;
		std::set<uintptr_t> edgeCareless;
		int curTI = Rcpp::as<Rcpp::IntegerVector>(travISet[curI])[fromQ];
		int nextI = curI + 1;
		while(nextI < maxTI){
			int nextTI = Rcpp::as<Rcpp::IntegerVector>(travISet[nextI])[fromQ];
			if(curTI != nextTI){ break; }
			nextI++;
		}
		GraphlineTraversalDump* toExP = Rcpp::as< Rcpp::XPtr<GraphlineTraversalDump> >(toExS[fromQ]).get();
		graphlineNoteMissedNodeEdge(sgrapP, &(toExP->travGraph), toExP->travEnds[curTI-1], &curHits, &edgeCareless);
		countTraversalsMissedNodesHelp(sgrapP, toExS, travISet, allCount, fromQ+1, curI, nextI, &curHits);
		curI = nextI;
	}
}

//' For each set of traversals, will count the number of missed nodes.
//'
//' @param sgrap The graph.
//' @param toExS The graph traversals.
//' @param travISet The sets of traversals.
//' @return The number of nodes missed by each set of traversals.
// [[Rcpp::export]]
Rcpp::IntegerVector countTraversalsMissedNodes(Rcpp::XPtr<GraphlineEJGraph> sgrap, Rcpp::List toExS, Rcpp::List travISet){
	std::vector<int> allCount;
	GraphlineEJGraph* sgrapP = sgrap.get();
	//check the traversals
	for(uintptr_t i = 0; i<toExS.size(); i++){
		GraphlineTraversalDump* toExP = Rcpp::as< Rcpp::XPtr<GraphlineTraversalDump> >(toExS[i]).get();
		if(toExP->graphAddr != (uintptr_t)sgrapP){ throw std::domain_error("Traversal did not come from the given graph."); }
	}
	for(uintptr_t i = 0; i<travISet.size(); i++){
		Rcpp::IntegerVector curTrav = Rcpp::as<Rcpp::IntegerVector>(travISet[i]);
		if(curTrav.size() != toExS.size()){ throw std::domain_error("Traversal indices do not match the traversal list."); }
		for(uintptr_t j = 0; j<toExS.size(); j++){
			int travInd = curTrav[j];
			GraphlineTraversalDump* toExP = Rcpp::as< Rcpp::XPtr<GraphlineTraversalDump> >(toExS[j]).get();
			if(travInd <= 0){ throw std::domain_error("Attempt to get index <= 0."); }
			if(travInd > toExP->travEnds.size()){ throw std::domain_error("Attempt to get traversal that does not exist."); }
		}
	}
	//run down the traversal enumeration
	std::set<uintptr_t> noNodes;
	countTraversalsMissedNodesHelp(sgrapP, toExS, travISet, &allCount, 0, 0, travISet.size(), &noNodes);
	//prepare the result
	return Rcpp::wrap(allCount);
}

/**
 * Run through all sets of individuals.
 * @param sgrapP The graph.
 * @param toExS The graph traversals for the individuals.
 * @param numInds The number of individuals in the mixture.
 * @param maxMissNode The maximum number of missed nodes to entertain.
 * @param curIndSet The current set of individuals.
 * @param allExpl The place to put all valid sets of individuals.
 * @param curI The current individual number (out of numInds).
 * @param minI The earliest individual to consider.
 */
void findExplainingIndividualsHelp(GraphlineEJGraph* sgrapP, Rcpp::List toExS, int numInds, int maxMissNode, std::vector<int>* curIndSet, std::vector< std::vector<int> >* allExpl, int curI, int minI){
	if(curI >= numInds){
		std::vector<int> curTIs;
			curTIs.insert(curTIs.end(), numInds, 0);
		while(true){
			std::set<uintptr_t> hitNodes;
			std::set<uintptr_t> hitEdges;
			for(int i = 0; i<numInds; i++){
				GraphlineTraversalDump* toExP = Rcpp::as< Rcpp::XPtr<GraphlineTraversalDump> >(toExS[(*curIndSet)[i]]).get();
				int curTI = curTIs[i];
				uintptr_t endTGI = toExP->travEnds[curTI];
				std::vector<GraphlineTraverseEntry>* curGrap = &(toExP->travGraph);
				std::set<uintptr_t>* hitNodeP = &hitNodes;
				std::set<uintptr_t>* hitEdgeP = &hitEdges;
				graphlineNoteMissedNodeEdge(sgrapP, curGrap, endTGI, hitNodeP, hitEdgeP);
			}
			if((sgrapP->allNodes.size() - hitNodes.size()) <= maxMissNode){
				allExpl->push_back(*curIndSet);
				break;
			}
			int fi = 0;
			while(fi < numInds){
				GraphlineTraversalDump* toExP = Rcpp::as< Rcpp::XPtr<GraphlineTraversalDump> >(toExS[(*curIndSet)[fi]]).get();
				curTIs[fi]++;
				if(curTIs[fi] < toExP->travEnds.size()){ break; }
				curTIs[fi] = 0;
				fi++;
			}
			if(fi >= numInds){ break; }
		}
		return;
	}
	for(int i = minI; i<toExS.size(); i++){
		if(Rcpp::as< Rcpp::XPtr<GraphlineTraversalDump> >(toExS[i]).get()->travEnds.size() == 0){ continue; }
		(*curIndSet)[curI] = i;
		findExplainingIndividualsHelp(sgrapP, toExS, numInds, maxMissNode, curIndSet, allExpl, curI+1, i+1);
	}
}


//' Find sets of individuals that can explain the graph.
//'
//' @param sgrap The graph.
//' @param toExS The graph traversals.
//' @param numInds The number of individuals to consider.
//' @param maxMissNode The maximum number of missed nodes.
//' @return The sets of individuals that can explain the graph.
// [[Rcpp::export]]
Rcpp::List findExplainingIndividuals(Rcpp::XPtr<GraphlineEJGraph> sgrap, Rcpp::List toExS, int numInds, int maxMissNode){
	GraphlineEJGraph* sgrapP = sgrap.get();
	//idiot checks
	if(numInds <= 0){ throw std::domain_error("Number of individuals must be positive."); }
	if(toExS.size() < numInds){
		Rcpp::List toRet;
		return toRet;
	}
	for(uintptr_t i = 0; i<toExS.size(); i++){
		GraphlineTraversalDump* toExP = Rcpp::as< Rcpp::XPtr<GraphlineTraversalDump> >(toExS[i]).get();
		if(toExP->graphAddr != (uintptr_t)sgrapP){ throw std::domain_error("Traversal did not come from the given graph."); }
	}
	//run down
	std::vector< std::vector<int> > allExpl;
	std::vector<int> curIndSet;
		curIndSet.insert(curIndSet.end(), numInds, 0);
	findExplainingIndividualsHelp(sgrapP, toExS, numInds, maxMissNode, &curIndSet, &allExpl, 0, 0);
	//report
	Rcpp::List toRet(allExpl.size());
	for(uintptr_t i = 0; i<allExpl.size(); i++){
		toRet[i] = Rcpp::wrap(allExpl[i]);
	}
	return toRet;
}

//' Get the edit sequence a traversal represents.
//'
//' @param sgrap The graph.
//' @param toEx The traversal set in question.
//' @param travInd The index of the traversal set.
//' @return The edit script to go from the query sequence to the reference graph. Three parallel vectors: the first is the location in the (current state of the) query, the second is the operation type (1 for substitution, 2 for insertion, 3 for deletion), and the third is the character (is garbage for deletion).
// [[Rcpp::export]]
Rcpp::List expandTraversalEditSequence(Rcpp::XPtr<GraphlineEJGraph> sgrap, Rcpp::XPtr<GraphlineTraversalDump> toEx, int travInd){
	GraphlineEJGraph* sgrapP = sgrap.get();
	GraphlineTraversalDump* toExP = toEx.get();
	if(travInd <= 0){ throw std::domain_error("Attempt to get index <= 0."); }
	if(travInd > toExP->travEnds.size()){ throw std::domain_error("Attempt to get traversal that does not exist."); }
	std::vector<int> changeLocs;
	std::vector<int> changeTypes;
	std::string changeChars;
	GraphlineTraverseEntry curEnt = toExP->travGraph[toExP->travEnds[travInd-1]];
	while(curEnt.prevEnt){
		GraphlineTraverseEntry prevEnt = toExP->travGraph[curEnt.prevEnt-1];
		if(curEnt.offCost != prevEnt.offCost){
			changeLocs.push_back(prevEnt.seqInd + 1);
			if(prevEnt.seqInd == curEnt.seqInd){
				changeTypes.push_back(2);
				changeChars.push_back(sgrap->allNodes[prevEnt.nodeInd].nodeChars[prevEnt.nodeSeqInd]);
			}
			else if(prevEnt.nodeSeqInd == curEnt.nodeSeqInd){
				changeTypes.push_back(3);
				changeChars.push_back(7);
			}
			else{
				changeTypes.push_back(1);
				changeChars.push_back(sgrap->allNodes[prevEnt.nodeInd].nodeChars[prevEnt.nodeSeqInd]);
			}
		}
		curEnt = prevEnt;
	}
	Rcpp::IntegerVector chanLocs = Rcpp::wrap(changeLocs);
	Rcpp::IntegerVector chanTypes = Rcpp::wrap(changeTypes);
	Rcpp::String chanChars = Rcpp::wrap(changeChars);
	Rcpp::List toRet(3);
	toRet[0] = chanLocs;
	toRet[1] = chanTypes;
	toRet[2] = chanChars;
	return toRet;
}

//' Get the variations from the reference that a traversal through the graph represents.
//'
//' @param refRS The reference sequence.
//' @param sgrap The graph.
//' @param toEx The traversal set in question.
//' @param travInd The index of the traversal set.
//' @return The variants for that traversal through the graph. Three parallel vectors: the first is the location in the reference, the second is the operation type (1 for substitution, 2 for insertion, 3 for deletion), and the third is the relevant character/sequence (is garbage for deletion).
// [[Rcpp::export]]
Rcpp::List expandTraversalGraphDifferences(Rcpp::String refRS, Rcpp::XPtr<GraphlineEJGraph> sgrap, Rcpp::XPtr<GraphlineTraversalDump> toEx, int travInd){
	//get all the nodes
	const char* refRSc = refRS.get_cstring();
	uintptr_t refRSLen = strlen(refRSc);
	GraphlineEJGraph* sgrapP = sgrap.get();
	GraphlineTraversalDump* toExP = toEx.get();
	if(travInd <= 0){ throw std::domain_error("Attempt to get index <= 0."); }
	if(travInd > toExP->travEnds.size()){ throw std::domain_error("Attempt to get traversal that does not exist."); }
	std::set<uintptr_t> allHitNode;
	GraphlineTraverseEntry curEnt = toExP->travGraph[toExP->travEnds[travInd-1]];
	allHitNode.insert(curEnt.nodeInd);
	while(curEnt.prevEnt){
		curEnt = toExP->travGraph[curEnt.prevEnt-1];
		allHitNode.insert(curEnt.nodeInd);
	}
	//get their variations
	std::vector<int> changeLocs;
	std::vector<int> changeTypes;
	std::vector<std::string> changeChars;
	for(std::set<uintptr_t>::iterator hitIt = allHitNode.begin(); hitIt != allHitNode.end(); hitIt++){
		GraphlineEJNode curNode = sgrapP->allNodes[*hitIt];
		std::pair<uintptr_t,uintptr_t> nodeRef = sgrapP->nodeRefCor[*hitIt];
		uintptr_t numRef = nodeRef.second - nodeRef.first;
		if(curNode.numChars == numRef){
			if(numRef){
				if(nodeRef.second > refRSLen){ throw std::domain_error("Reference provided was not used in the creation of the graph."); }
				if(memcmp(refRSc + nodeRef.first, curNode.nodeChars, numRef) != 0){
					changeLocs.push_back(nodeRef.first);
					changeTypes.push_back(1);
					changeChars.push_back(std::string(curNode.nodeChars, curNode.nodeChars + curNode.numChars));
				}
			}
		}
		else if(curNode.numChars < numRef){
			changeLocs.push_back(nodeRef.first);
			changeTypes.push_back(3);
			changeChars.push_back((std::string()));
		}
		else if(curNode.numChars > numRef){
			changeLocs.push_back(nodeRef.first);
			changeTypes.push_back(2);
			changeChars.push_back(std::string(curNode.nodeChars, curNode.nodeChars + curNode.numChars));
		}
	}
	Rcpp::IntegerVector chanLocs = Rcpp::wrap(changeLocs);
	Rcpp::IntegerVector chanTypes = Rcpp::wrap(changeTypes);
	Rcpp::StringVector chanChars = Rcpp::wrap(changeChars);
	Rcpp::List toRet(3);
	toRet[0] = chanLocs;
	toRet[1] = chanTypes;
	toRet[2] = chanChars;
	return toRet;
}
