// For compilers that support precompilation, includes "wx/wx.h".
#include "wx/wxprec.h"
#include <wx/arrimpl.cpp>


#ifdef __BORLANDC__
#pragma hdrstop
#endif

#ifndef WX_PRECOMP
#include <wx/wx.h>
#endif

#include <wx/image.h>
#include <wx/file.h>
#include <wx/bitmap.h>
#include <wx/textfile.h>

#include <vector>
#include <iostream>
#include <fstream>
#include <queue>
#include <deque>
#include <cmath>
#include <memory>
#include <map>
#include <chrono>
#include <boost/histogram.hpp>
#include <boost/format.hpp>
using namespace boost::histogram;
using namespace std;
using namespace std::chrono;

using std::vector;
using std::cout;
using std::endl;
using std::queue;
using std::deque;
using std::min;
using std::ifstream;
#define APP_NAME "Image Segmentation"

unsigned int _SIGMA = 2;
unsigned int _LAMBDA = 0;

const bool _METRICS = false; // true for testing metrics

bool maxFlowInitalComputed = false;
//unsigned int bkgTimes = 0;
//unsigned int objTimes = 0;

vector <int> bkgDots;
vector <int> objDots;
int maxProb = 0;

double firstMetric(unsigned char* frame1, unsigned char* frame2)
{
    unsigned int correct = 0;
    int w = 608;
    int h = 511;

    int intensity1, intensity2;

    for (int i = 0; i < h; i++)
    {
        for (int j = 0; j < w; j++)
        {
            // for a single channel grey scale image
            intensity1 = frame1[i * w + j];
            intensity2 = frame2[i * w + j];
            if (intensity1 == intensity2) correct += 1;
        }
    }
    //cout << "correct: " << correct << endl;
    //cout << "all: " << frame1->cols * frame1->rows << endl;
    return correct / double(w*h);
}

double jaccardMetric(unsigned char* frame1, unsigned char* frame2) // for object 
{
    unsigned int intersect = 0;
    unsigned int unionsect = 0;

    int w = 608;
    int h = 511;

    int intensity1, intensity2;

    for (int i = 0; i < h; i++)
    {
        for (int j = 0; j < w; j++)
        {
            // for a single channel grey scale image
            intensity1 = frame1[i * w + j];
            intensity2 = frame2[i * w + j];
            if (intensity1 == intensity2 && intensity2 == 255) intersect += 1;
            if (intensity1 == 255 || intensity2 == 255) unionsect += 1;
        }
    }

    //cout << "\nintersect: " << intersect << endl;
    //cout << "unionsect: " << unionsect << endl;
    //cout << "all:       " << frame1->cols * frame1->rows << endl;
    return intersect / double(unionsect);
}

struct Histogram
{
    Histogram(size_t n, int low, double high)
    {
        int x = low;
        double dx = (high - low) / n;
        while (x <= high) { _bins[x] = 0; x += dx; }
    }
    void fill(std::vector <int> &values)
    {
        for (int val : values)
        {
            if (val < _bins.begin()->first || val > _bins.rbegin()->first)
                return;

            std::prev(_bins.upper_bound(val))->second++;
        }
    }

    int at(int depth) // return the number of values in bin for depth=value
    {
        if (depth == 0) return _bins[0];
        return (std::prev(_bins.lower_bound(depth)))->second;
    }

    void print(std::ostream& o) const
    {
        for (auto b1 = _bins.begin(), b2 = std::next(b1);
            b2 != _bins.end(); ++b1, ++b2)
            o << b1->first << "\t- " << b2->first << ":\t"
            << std::string(b1->second, '*')
            << std::endl;
    }
    std::map<int, size_t> _bins;
};

struct Edge
{
    // To store current flow and capacity of edge
    double flow, capacity;

    // An edge u--->v has start vertex as u and end
    // vertex as v.
    int u, v;

    int pair; // adj[v][pair] == adj[vu]

    Edge(double flow, double capacity, int u, int v) : flow(flow), capacity(capacity), u(u), v(v), pair(-1) {}
};

// Represent a Vertex
struct Vertex
{
    int name, h, depth, pixelX, pixelY, label;
    double e_flow;

    Vertex(int name, int h, double e_flow, 
        int x, int y, int label = 2, int depth = -1) :
       name(name),  h(h), e_flow(e_flow), depth(depth), pixelX(x), pixelY(y), label(label) {}
};



double bValue(Vertex* p, Vertex* q, int sigma = _SIGMA, int dist = 1)
{
    return exp(-pow(p->depth - q->depth, 2) / (2 * pow(sigma, 2))) / dist;
}


// To represent a flow network
class Graph
{
    int E; // No. of edges 

    vector<Vertex> ver;
    // adjacency list of residual graph
    vector< vector<Edge> > adj;

    queue <Vertex*> active;

    double probabilityValueMax;
    double maxB;

    // Function to push excess flow from u
    bool push(int u);

    // Function to relabel a vertex u
    void relabel(int u);

    // This function is called to initialize
    // preprocess
    void preprocess(int s);

    // Function to reverse edge
    void updateReverseEdgeFlow(int uAdj, int i, double flow);

public:

    int V; // No. of vertices

    Graph(int V) : objHist(nullptr), bkgHist(nullptr), probabilityValueMax(0), maxB(0.0){
        this->V = V;
        this->E = 0;

        // all vertices are initialized with 0 height
        // and 0 excess flow
        for (int i = 0; i < V; i++) {
            ver.push_back(Vertex(i, 0, 0, 0, 0, 2));
            adj.push_back(vector<Edge>());
        }
    }

    void setLabelVertex(int positionInVer, int label)
    {
        ver[positionInVer].label = label;
    }

    int getV() { return V; }

    int getLabel(int verPosition)
    {
        return ver[verPosition].label;
    }
    
    // function to add an edge to graph
    void addEdge(int u, int v, double capacity)
    {
        // flow is initialized with 0 for all edge

        int index = findEdge(u, v);
        if (index == -1) // there is no u<->v in adj list 
        {
            adj[u].push_back(Edge(capacity, capacity, u, v));
            adj[v].push_back(Edge(capacity, capacity, v, u));

            (adj[u][adj[u].size() - 1]).pair = adj[v].size() - 1;
            (adj[v][adj[v].size() - 1]).pair = adj[u].size() - 1;

            E += 1; // as in non-directed graph

        }
    }

    void updateEdge(int u, int area)
    {
        int lastPos = adj[u + 1].size() - 1;
        double cP = adj[u + 1][lastPos-1].flow + adj[u + 1][lastPos].flow; // 1 - with Source, 2 - with Sink

        if (area) // bkg
        {
            adj[u + 1][lastPos - 1].capacity = cP; // with source
            adj[0][adj[u + 1][lastPos - 1].pair].capacity = cP;

            adj[u + 1][lastPos].capacity = 1 + maxB + cP; // with sink
            adj[V - 1][adj[u + 1][lastPos].pair].capacity = 1 + maxB + cP;;
        }
        else // obj
        {
            adj[u + 1][lastPos - 1].capacity = 1 + maxB + cP; // with source
            adj[0][adj[u + 1][lastPos - 1].pair].capacity = 1 + maxB + cP;;

            adj[u + 1][lastPos].capacity = cP; // with sink
            adj[V - 1][adj[u + 1][lastPos].pair].capacity = cP;
        }
    }

    // Return number in adj[u] for edge u->v
    int findEdge(int u, int v)
    {
        for (unsigned int i = 0; i < adj[u].size(); ++i)
        {
            if (adj[u][i].v == v)
            {
                return i;
            }
        }
        return -1;
    }

    void globalRelabeling();
    void bfs(int start, vector<int>& distances);

    // returns maximum flow from s to t
    double getMaxFlow(int s);
    void minCut(vector <Vertex*>& mincut);

    void Init(int numColumn, int numRow, unsigned char* imageArray);
    void contactWithTerminals();

    std::unique_ptr<Histogram> objHist;
    std::unique_ptr<Histogram> bkgHist;

    void computeHistogram()
    {
        objHist = std::make_unique<Histogram>(log(objDots.size()) + 1, 0, 255);
        bkgHist = std::make_unique<Histogram>(log(bkgDots.size()) + 1, 0, 255);
        objHist->fill(objDots);
        bkgHist->fill(bkgDots);
    }

    double probabilityValue(int depth, int area, int lambda=_LAMBDA);

    void PrintCondition()
    {
        for (unsigned int i = 0; i < ver.size(); ++i)
            cout << "Node " << ver[i].name << " with excess " << ver[i].e_flow << " and distance " << ver[i].h << endl;
        for (unsigned int i = 0; i < adj.size(); ++i)
        {
            for (unsigned int j = 0; j < adj[i].size(); ++j)
                cout << "Edge " << adj[i][j].u << "-->" << adj[i][j].v << ", " << adj[i][j].flow << "/"
                << adj[i][j].capacity << "\n";
        }
        cout << "\n\n";
    }
};

void Graph::Init(int numColumn, int numRow, unsigned char* imageArray)
{
    /* INMPORTANT : ver[0] == source, then image's vertices starting from the 1 index */
    int numVer = numColumn * numRow;
    int currentPixelPosition = 0; int verPosition = 0;
    int depthPixel = 0; // r=g=b
    

    // WASD - neighbours pixels (in ver-vector)
    int nA = 0;
    int nW = 0;
    int nS = 0;
    int nD = 0;

    vector <double> bValues;

    for (int i = 0; i < numRow; ++i)
    {
        for (int j = 0; j < numColumn; ++j)
        {
            currentPixelPosition = numColumn * i + j; // in image => m_myImage[3*curPP] 
            verPosition = currentPixelPosition + 1; // +1 due to [source = 0 index and name]
            nA = currentPixelPosition-1; nW = currentPixelPosition - numColumn; 
            nS = currentPixelPosition + numColumn; nD = currentPixelPosition + 1;

            depthPixel = imageArray[3 * currentPixelPosition];
            ver[verPosition].depth = depthPixel;
            ver[verPosition].pixelX = i;
            ver[verPosition].pixelY = j;

            if (i == 0)
            {
                if (j == 0)                                     // г nS nD
                {
                    ver[nS+1].depth = imageArray[3 * nS];
                    ver[nD+1].depth = imageArray[3 * nD];
                    bValues.push_back(bValue(&ver[verPosition], &ver[nS + 1]));
                    bValues.push_back(bValue(&ver[verPosition], &ver[nD + 1]));
                    maxB = maxB < *max_element(begin(bValues), end(bValues)) ? *max_element(begin(bValues), end(bValues)) : maxB;
                    addEdge(verPosition, nS + 1, bValues[0]);
                    addEdge(verPosition, nD + 1, bValues[1]);
                    bValues.clear();
                    continue; }
                else if (j == numColumn - 1)                    // ˥ nA nS
                {
                    ver[nA + 1].depth = imageArray[3 * nA];
                    ver[nS + 1].depth = imageArray[3 * nS];
                    bValues.push_back(bValue(&ver[verPosition], &ver[nA + 1]));
                    bValues.push_back(bValue(&ver[verPosition], &ver[nS + 1]));
                    maxB = maxB < *max_element(begin(bValues), end(bValues)) ? *max_element(begin(bValues), end(bValues)) : maxB;
                    addEdge(verPosition, nA + 1, bValues[0]);
                    addEdge(verPosition, nS + 1, bValues[1]);
                    bValues.clear();
                    continue; }
                else                                            // T nA nS nD 
                {
                    ver[nA + 1].depth = imageArray[3 * nA];
                    ver[nD + 1].depth = imageArray[3 * nD];
                    ver[nS + 1].depth = imageArray[3 * nS];
                    bValues.push_back(bValue(&ver[verPosition], &ver[nA + 1]));
                    bValues.push_back(bValue(&ver[verPosition], &ver[nD + 1]));
                    bValues.push_back(bValue(&ver[verPosition], &ver[nS + 1]));
                    maxB = maxB < *max_element(begin(bValues), end(bValues)) ? *max_element(begin(bValues), end(bValues)) : maxB;
                    addEdge(verPosition, nA + 1, bValues[0]);
                    addEdge(verPosition, nD + 1, bValues[1]);
                    addEdge(verPosition, nS + 1, bValues[2]);
                    bValues.clear();
                    continue; }
            }
            else if (i == numRow - 1)
            {
                if (j == 0)                                     // ˪ nW nD
                {
                    ver[nD + 1].depth = imageArray[3 * nD];
                    ver[nW + 1].depth = imageArray[3 * nW];
                    bValues.push_back(bValue(&ver[verPosition], &ver[nD + 1]));
                    bValues.push_back(bValue(&ver[verPosition], &ver[nW + 1]));
                    maxB = maxB < *max_element(begin(bValues), end(bValues)) ? *max_element(begin(bValues), end(bValues)) : maxB;
                    addEdge(verPosition, nD + 1, bValues[0]);
                    addEdge(verPosition, nW + 1, bValues[1]);
                    bValues.clear();
                    continue; }
                else if (j == numColumn - 1)                    // ˩ nA nW
                {
                    ver[nA + 1].depth = imageArray[3 * nA];
                    ver[nW + 1].depth = imageArray[3 * nW];
                    bValues.push_back(bValue(&ver[verPosition], &ver[nA + 1]));
                    bValues.push_back(bValue(&ver[verPosition], &ver[nW + 1]));
                    maxB = maxB < *max_element(begin(bValues), end(bValues)) ? *max_element(begin(bValues), end(bValues)) : maxB;
                    addEdge(verPosition, nA + 1, bValues[0]);
                    addEdge(verPosition, nW + 1, bValues[1]);
                    bValues.clear();
                    continue; }
                else                                            // ⊥ nA nW nD
                {
                    ver[nA + 1].depth = imageArray[3 * nA];
                    ver[nD + 1].depth = imageArray[3 * nD];
                    ver[nW + 1].depth = imageArray[3 * nW];
                    bValues.push_back(bValue(&ver[verPosition], &ver[nA + 1]));
                    bValues.push_back(bValue(&ver[verPosition], &ver[nD + 1]));
                    bValues.push_back(bValue(&ver[verPosition], &ver[nW + 1]));
                    maxB = maxB < *max_element(begin(bValues), end(bValues)) ? *max_element(begin(bValues), end(bValues)) : maxB;
                    addEdge(verPosition, nA + 1, bValues[0]);
                    addEdge(verPosition, nD + 1, bValues[1]);
                    addEdge(verPosition, nW + 1, bValues[2]);
                    bValues.clear();
                    continue; }
            }
            else if (j == 0) {                                  // left border nW nD nS
                ver[nS+ 1].depth = imageArray[3 * nS];
                ver[nD + 1].depth = imageArray[3 * nD];
                ver[nW + 1].depth = imageArray[3 * nW];
                bValues.push_back(bValue(&ver[verPosition], &ver[nS + 1]));
                bValues.push_back(bValue(&ver[verPosition], &ver[nD + 1]));
                bValues.push_back(bValue(&ver[verPosition], &ver[nW + 1]));
                maxB = maxB < *max_element(begin(bValues), end(bValues)) ? *max_element(begin(bValues), end(bValues)) : maxB;
                addEdge(verPosition, nS + 1, bValues[0]);
                addEdge(verPosition, nD + 1, bValues[1]);
                addEdge(verPosition, nW + 1, bValues[2]);
                bValues.clear();
                continue; } 
            else if (j == numColumn - 1) {                      // ˧ nW nA nS
                ver[nA + 1].depth = imageArray[3 * nA];
                ver[nW + 1].depth = imageArray[3 * nW];
                ver[nS + 1].depth = imageArray[3 * nS];
                bValues.push_back(bValue(&ver[verPosition], &ver[nA + 1]));
                bValues.push_back(bValue(&ver[verPosition], &ver[nW + 1]));
                bValues.push_back(bValue(&ver[verPosition], &ver[nS + 1]));
                maxB = maxB < *max_element(begin(bValues), end(bValues)) ? *max_element(begin(bValues), end(bValues)) : maxB;
                addEdge(verPosition, nA + 1, bValues[0]);
                addEdge(verPosition, nW + 1, bValues[1]);
                addEdge(verPosition, nS + 1, bValues[2]);
                bValues.clear();
                continue; }      
            else {                                              // + nW nA nS nD
                ver[nA + 1].depth = imageArray[3 * nA];
                ver[nW + 1].depth = imageArray[3 * nW];
                ver[nS + 1].depth = imageArray[3 * nS];
                ver[nD + 1].depth = imageArray[3 * nD];
                bValues.push_back(bValue(&ver[verPosition], &ver[nA + 1]));
                bValues.push_back(bValue(&ver[verPosition], &ver[nW + 1]));
                bValues.push_back(bValue(&ver[verPosition], &ver[nS + 1]));
                bValues.push_back(bValue(&ver[verPosition], &ver[nD + 1]));
                maxB = maxB < *max_element(begin(bValues), end(bValues)) ? *max_element(begin(bValues), end(bValues)) : maxB;
                addEdge(verPosition, nA + 1, bValues[0]);
                addEdge(verPosition, nW + 1, bValues[1]);
                addEdge(verPosition, nS + 1, bValues[2]);
                addEdge(verPosition, nD + 1, bValues[3]);
                bValues.clear();
                continue; }                              
        }
    }

    // source and sink has no edges with image's graph for now
}

void Graph::contactWithTerminals()
{
    // arcs beetwen pixel-vertex and source/sink
    int source = 0;
    int sink = V - 1;

    for (int i = 1; i < V - 1; ++i)
    {
        if (ver[i].label == 0) // obj dot
        {
            addEdge(i, source, 1 + maxB);
            addEdge(i, sink, 0);
        }
        else if (ver[i].label == 1) // bkg dot
        {
            addEdge(i, source, 0);
            addEdge(i, sink, 1 + maxB); // 1 + ver[i].maxB == K
        }
        else // unknown dot
        {
            addEdge(i, source, probabilityValue(ver[i].depth, 1)); // bkg
            addEdge(i, sink, probabilityValue(ver[i].depth, 0));   // obj
        }
    }
}

double Graph::probabilityValue(int depth, int area, int lambda) // area: 0 - obj, bkg - 1; int lambda = 50
{
    double probObj = double(objHist->at(depth)) / objDots.size();
    double probBkg = double(bkgHist->at(depth)) / bkgDots.size();
    double sumProb = probObj + probBkg;

    probabilityValueMax = probObj > probabilityValueMax ? probObj : probabilityValueMax;
    probabilityValueMax = probBkg > probabilityValueMax ? probBkg : probabilityValueMax;

    if (sumProb < 1.0e-10 || probObj < 1.0e-10 || probBkg < 1.0e-10) //==0
    {
        return 100*probabilityValueMax; // as double max
    }

    // with normalization
    if (area) // bkg
    {
        return -lambda * log(probBkg / sumProb);
    }
    else // obj
    {
        return -lambda * log(probObj / sumProb);
    }
}

void Graph::preprocess(int s)
{
    // Making h of source Vertex equal to no. of vertices
    // Height of other vertices is 0.
    ver[s].h = ver.size();
    double flow;
    Vertex endV(0, 0, 0, 0, 0);

    //
    for (unsigned int i = 0; i < adj[s].size(); i++)
    {
        endV = ver[adj[s][i].v];

        flow = adj[s][i].capacity;
        adj[s][i].flow = 0;
        adj[s][i].capacity = 0; // means not having this edge in residual

        // Initialize excess flow for adjacent v
        ver[endV.name].e_flow += flow;
        if (ver[endV.name].name != ver.size() - 1)
        {
            active.push(&ver[endV.name]);
        }

        // update for s<-v (+flow to capacity)
        updateReverseEdgeFlow(s, i, flow);
    }
}

// Update reverse flow for flow added on ith Edge for u
void Graph::updateReverseEdgeFlow(int uAdj, int i, double flow)
{
    //    cout << "Update reverse for " << uAdj << "-->" << adj[uAdj][i].v << "with flow " << flow << endl;
    int u = adj[uAdj][i].v;
    adj[u][adj[uAdj][i].pair].capacity += flow;
}

// To push flow from overflowing vertex u
bool Graph::push(int u)
{
        // Traverse through all edges to find an adjacent (of u)
        // to which flow can be pushed
    for (unsigned int i = 0; i < adj[u].size(); i++)
    {
        // if capacity==0 then no flow can be pushed
        if (adj[u][i].capacity == 0)
            continue;

        // Push is only possible if height of adjacent
        // is smaller than height of overflowing vertex
        if (ver[u].h > ver[adj[u][i].v].h && ver[u].e_flow > 0)
        {

            // Flow to be pushed is equal to minimum of
            // residual capacity of edge and excess flow.
            double flow = min(adj[u][i].capacity,
                ver[u].e_flow);


            // Reduce excess flow for overflowing vertex
            ver[u].e_flow -= flow;

            // Increase excess flow for adjacent
            ver[adj[u][i].v].e_flow += flow;


            // If v!= source and !=sink
            if (adj[u][i].v && adj[u][i].v != ver.back().name)
            {
                // and if excess in v was 0, then add in active
                if (!ver[adj[u][i].v].e_flow - flow)
                {
                   active.push(&ver[adj[u][i].v]);
                }
            }


            // Add residual flow (With capacity 0 and negative flow)
            adj[u][i].capacity -= flow;

            updateReverseEdgeFlow(u, i, flow);
            return true;
        }
    }
    return false;
}

// function to relabel vertex u
void Graph::relabel(int u)
{
    //    cout << "\n RELABEL\n";
        // Initialize minimum height of an adjacent
    int mh = INT_MAX;

    // Pop the vertex from the queue

//    cout << "Pop from \"active\" node " << active.front()->name << " with excess " << active.front()->e_flow << endl;
    active.pop();

    // Find the adjacent with minimum height
    for (unsigned int i = 0; i < adj[u].size(); i++)
    {
        // if flow is equal to capacity then no
        // relabeling
        if (adj[u][i].capacity == 0)
            continue;

        // Update minimum height
        if (ver[adj[u][i].v].h < mh)
        {
            mh = ver[adj[u][i].v].h;

            // updating height of u
            ver[u].h = mh + 1;
        }
    }

    // If excess if equal to 0, then no pushing in active
    if (!ver[u].e_flow) return;

    if (ver[u].name != ver.size() - 1) // СЃС‚РѕРє РЅРµ СЃС‡РёС‚Р°РµС‚СЃСЏ Р°РєС‚РёРІРЅРѕР№ РІРµСЂС€РёРЅРѕР№
    {
        active.push(&ver[u]);
        //        cout << "Push in \"active\" node " << ver[u].name << " with excess " << active.back()->e_flow << endl;
    }
    //    PrintCondition();
}

void Graph::minCut(vector <Vertex*>& mincut)
{
    int start = 0;
    vector <bool> visited(V, 0);

    deque <int> deq;
    deq.push_back(start);
    mincut.push_back(&ver[start]);
    visited[start] = true;

    int neighbour;
    while (!deq.empty())
    {
        int u = deq.front();
        deq.pop_front();
        for (int i = 0; i < adj[u].size(); ++i)
        {
            neighbour = adj[u][i].v;
            //            cout << "adj[" << u << "][" << i << "].capacity = " << fixed << adj[u][i].capacity << endl;
            if (adj[u][i].capacity > 0.0 && visited[neighbour] == false)
            {
                //                cout << " and its unvisited neighbour : " << neighbour << endl;
                deq.push_back(neighbour);
                mincut.push_back(&ver[neighbour]);
                ver[neighbour].label = 0;
                visited[neighbour] = true;
            }
        }
    }

}


void Graph::bfs(int start, vector<int>& distances)
{
    //distances[start] = 0;
    vector <bool> visited(V, 0);

    // not queue because of the easier printing the elements
    deque <int> deq;
    deq.push_back(start);

    visited[start] = true;
    // source is 0 vertex
    distances[start] = start == 0 ? ver[0].h : 0;

    int coutG = 0;
    int neighbour;
    while (!deq.empty())
    {
        int u = deq.front();
        deq.pop_front();

        for (int i = 0; i < adj[u].size(); ++i)
        {
            neighbour = adj[u][i].v;

            // found residual for u->(v) through v : adj[neighbour][adj[u][i].pair] = adj[vu]
            if (adj[neighbour][adj[u][i].pair].capacity != 0 && visited[neighbour] == false)
            {
                distances[neighbour] = distances[u] + 1;
                deq.push_back(neighbour);
                visited[neighbour] = true;
            }
        }
    }
}

void Graph::globalRelabeling()
{
    /*
        Запускаем бфс от каждой вершины сети, находим расстояние для неё до стока по остаточным ребрам, меняем высоту.
        Если же добраться до стока невозможно, то находим расстояние от истока до вершины и заменяем высоту на него.
     */

    int terminal = ver.size() - 1;
    int source = 0;

    vector <int> sourceD(V, 0);
    vector <int> sinkD(V, 0);

    bfs(terminal, sinkD);
    bfs(source, sourceD); // до какой вершины (source|sink), массив для заполнения расстояниями
    for (int i = 1; i < terminal; ++i)
    {
        ver[i].h = sinkD[i];
        if (!ver[i].h) ver[i].h = sourceD[i];
    }

}


// main function for printing maximum flow of graph
double Graph::getMaxFlow(int s)
{
    int countGB = 0;
    preprocess(s);

    ofstream log;
    log.open("log.txt", std::ios_base::app);


    // loop until none of the Vertex is active
    while (!active.empty())
    {
        int u = active.front()->name;
        if (!push(u))
        {
            if (countGB % (V + E)) //countGB%(3*V)
            {
                relabel(u);
                countGB += 1;
            }
            else
            {
                //                cout << "else Gb=" << countGB << endl;


                auto start = chrono::high_resolution_clock::now();

                globalRelabeling();

                auto end = chrono::high_resolution_clock::now();
                double time_taken = chrono::duration_cast<chrono::nanoseconds>(end - start).count();
                time_taken *= 1e-9;
                log << time_taken << "\n\n";

                countGB += 1;
            }
        }
    }

    log.close();
    // ver.back() returns last Vertex, whose
    // e_flow will be final maximum flow
    return ver.back().e_flow;
}


enum
{
	ID_QUIT = 1,
	ID_ABOUT,
	ID_LOAD,
	ID_SAVE,
	ID_PROCESS,
	ID_BEST_SIZE
};

//************************************************************************
//************************************************************************
// Canvas class (where we display the image)
//************************************************************************
//************************************************************************

//------------------------------------------------------------------------
class MyCanvas : public wxPanel
	//------------------------------------------------------------------------
{
public:
	MyCanvas(wxWindow* parent, wxWindowID, const wxPoint& pos, const wxSize& size);
	~MyCanvas();

	void LoadImage(wxString fileName);
	void SaveImage(wxString fileName);
	void ProcessImage();
	void BestSize();

	wxStaticText* st1;
	wxStaticText* st2;
    wxStaticText* st3;

    void OnLeftMove(wxMouseEvent& event);
    void OnLeftDown(wxMouseEvent& event); // mouse event handler
    void OnLeftUp(wxMouseEvent& event);    
    
    void OnRightDown(wxMouseEvent& event);
    void OnRightUp(wxMouseEvent& event);

protected:
    bool isLeftDown;
    bool isRightDown;
    wxPoint penPos;


private:
	int m_imageWidth;
	int m_imageHeight;
	wxBitmap m_imageBitmap;	// used to display the image
	wxImage* m_imageRGB;		// used to load the image
	unsigned char* m_myImage;	// used to process the image

    std::unique_ptr<Graph> p_V; // pointer on representation of image via graph

	void OnPaint(wxPaintEvent& event);

    friend class Graph;

	DECLARE_EVENT_TABLE()
};

BEGIN_EVENT_TABLE(MyCanvas, wxPanel)
    EVT_PAINT(MyCanvas::OnPaint)
    EVT_LEFT_DOWN(MyCanvas::OnLeftDown)
    EVT_LEFT_UP(MyCanvas::OnLeftUp)
    EVT_RIGHT_DOWN(MyCanvas::OnRightDown)
    EVT_RIGHT_UP(MyCanvas::OnRightUp)
    EVT_MOTION(MyCanvas::OnLeftMove)
END_EVENT_TABLE()



//------------------------------------------------------------------------
MyCanvas::MyCanvas(wxWindow* parent, wxWindowID id,
	const wxPoint& pos, const wxSize& size)
	: wxPanel(parent, id, pos, size, wxSUNKEN_BORDER)
	//------------------------------------------------------------------------
{
	st1 = new wxStaticText(this, -1, wxT("X"), wxPoint(10, 10));
	st2 = new wxStaticText(this, -1, wxT("Y"), wxPoint(10, 30));
    st3 = new wxStaticText(this, -1, wxT("Z"), wxPoint(10, 50));
	m_myImage = NULL;
	m_imageRGB = NULL;
    isLeftDown = false;
    isRightDown = false;

}

//------------------------------------------------------------------------
MyCanvas::~MyCanvas()
//------------------------------------------------------------------------
{
	if (m_myImage)
		free(m_myImage);
	if (m_imageRGB)
		delete m_imageRGB;
}

//------------------------------------------------------------------------
void MyCanvas::LoadImage(wxString fileName)
//------------------------------------------------------------------------
{
	if (m_myImage)
		free(m_myImage);
	if (m_imageRGB)
		delete m_imageRGB;

    //if (p_V) p_V.reset();

	// open image dialog box
	m_imageRGB = new wxImage(fileName, wxBITMAP_TYPE_ANY, -1); // ANY => can load many image formats
	m_imageBitmap = wxBitmap(*m_imageRGB, -1); // ...to get the corresponding bitmap

	m_imageWidth = m_imageRGB->GetWidth();
	m_imageHeight = m_imageRGB->GetHeight();

	//m_myImage = (unsigned char*)malloc(m_imageWidth * m_imageHeight * 3);
	//memcpy(m_myImage, m_imageRGB->GetData(), m_imageWidth * m_imageHeight * 3);
    //st1->SetLabel(wxString::Format(wxT("m_imageWidth = %d"), m_imageWidth));
    //st2->SetLabel(wxString::Format(wxT("m_imageHeight = %d"), m_imageHeight));


    /*
    Initialization of graph.
    */
    //p_V = std::make_unique<Graph>(m_imageWidth* m_imageHeight + 2); //+2 for sink and source

    if (!_METRICS)
    {
        m_myImage = (unsigned char*)malloc(m_imageWidth * m_imageHeight * 3);
        memcpy(m_myImage, m_imageRGB->GetData(), m_imageWidth * m_imageHeight * 3);
        st1->SetLabel(wxString::Format(wxT("m_imageWidth = %d"), m_imageWidth));
        st2->SetLabel(wxString::Format(wxT("m_imageHeight = %d"), m_imageHeight));

        p_V = std::make_unique<Graph>(m_imageWidth * m_imageHeight + 2); //+2 for sink and source
        p_V->Init(m_imageWidth, m_imageHeight, m_myImage);
        st2->SetLabel(wxString::Format(wxT("V = %d"), p_V->getV()));

    }
    else
    {
        ofstream log;
        log.open("param_test1.txt", std::ios_base::app);

        vector<unsigned int> sigmas = {110, 120, 140, 250}; //  1, 2, 3, 4, 10, 20, 40, 50, 60, 80, 100    };
        vector <unsigned int> lambdas = { 0, 1, 2, 3, 4, 5, 10, 20, 30, 40, 50, 60 };

        for (unsigned int sigma : sigmas)
        {
            for (unsigned int lambda : lambdas)
            {
                _SIGMA = sigma;
                _LAMBDA = lambda;

                m_myImage = (unsigned char*)malloc(m_imageWidth * m_imageHeight * 3);
                memcpy(m_myImage, m_imageRGB->GetData(), m_imageWidth * m_imageHeight * 3);
                st1->SetLabel(wxString::Format(wxT("sigma = %d"), sigma));
                st2->SetLabel(wxString::Format(wxT("lambda = %d"), lambda));


                p_V = std::make_unique<Graph>(m_imageWidth * m_imageHeight + 2);
                p_V->Init(m_imageWidth, m_imageHeight, m_myImage);
                //st2->SetLabel(wxString::Format(wxT("V = %d"), p_V->getV()));

                int n, m, x, y, c, verPos;
                ifstream bkgfile("C:\\Users\\Asus\\CLionProjects\\study\\graphProject\\bkg_test.txt");
                bkgfile >> n >> m;
                while (bkgfile >> x >> y >> c)
                {
                    verPos = y * n + x;
                    p_V->setLabelVertex(verPos + 1, 1);
                    bkgDots.push_back(m_myImage[3 * verPos]);
                }

                ifstream objfile("C:\\Users\\Asus\\CLionProjects\\study\\graphProject\\obj_test.txt");
                objfile >> n >> m;
                while (objfile >> x >> y >> c)
                {
                    verPos = y * n + x;
                    p_V->setLabelVertex(verPos + 1, 0);
                    objDots.push_back(m_myImage[3 * verPos]);
                }

                p_V->computeHistogram();
                p_V->contactWithTerminals();

                auto start = chrono::high_resolution_clock::now();

                st1->SetLabel(wxString::Format(wxT("MaxFlow = %f"), p_V->getMaxFlow(0)));
                maxFlowInitalComputed = true;

                auto end = chrono::high_resolution_clock::now();
                double time_taken =
                    chrono::duration_cast<chrono::nanoseconds>(end - start).count();

                time_taken *= 1e-9;
                //st2->SetLabel(wxString::Format(wxT("Time MF: %f"), time_taken));

                vector <Vertex*> mincut;
                p_V->minCut(mincut);
                st3->SetLabel(wxString::Format(wxT("mincut : %d/%d"), mincut.size(), p_V->getV()));


                int label = 2;
                int position = 0;

                for (unsigned int j = 0; j < m_imageHeight; ++j)
                {
                    for (unsigned int k = 0; k < m_imageWidth; ++k)
                    {
                        position = m_imageWidth * j + k;
                        label = p_V->getLabel(position + 1);
                        m_myImage[3 * position] = label == 0 ? 255 : 0;
                        m_myImage[3 * position + 1] = label == 0 ? 255 : 0;
                        m_myImage[3 * position + 2] = label == 0 ? 255 : 0;
                    }
                }


                wxImage* ideal_imageRGB = new wxImage("C:\\Users\\Asus\\CLionProjects\\study\\graphProject\\test_ideal.bmp", wxBITMAP_TYPE_ANY, -1);
                unsigned char* ideal_myImage = (unsigned char*)malloc(m_imageWidth * m_imageHeight * 3);
                memcpy(ideal_myImage, ideal_imageRGB->GetData(), ideal_imageRGB->GetWidth() * ideal_imageRGB->GetHeight() * 3);

                log << _SIGMA << " " << _LAMBDA << " " << firstMetric(m_myImage, ideal_myImage) << " " << jaccardMetric(m_myImage, ideal_myImage) << "\n";

            }
        }

        log.close();
    }

    

	// update GUI size
	SetSize(m_imageWidth, m_imageHeight);
	GetParent()->SetClientSize(GetSize());

	// update display
	Refresh(false);
}

//------------------------------------------------------------------------
void MyCanvas::SaveImage(wxString fileName)
//------------------------------------------------------------------------
{
	bool b;

	wxImage* tempImage = new wxImage(m_imageWidth, m_imageHeight, m_myImage, true); // lend my image buffer...
	b = tempImage->SaveFile(fileName);
	delete(tempImage);		// buffer not needed any more

	if (!b)
		wxMessageBox(wxT("A problem occured during saving"));
}

//------------------------------------------------------------------------
void MyCanvas::ProcessImage()
{
	long int i = m_imageWidth * m_imageHeight * 3;

    if (!maxFlowInitalComputed)
    {
        p_V->computeHistogram();
        p_V->contactWithTerminals();
    }


auto start = chrono::high_resolution_clock::now();

    st1->SetLabel(wxString::Format(wxT("MaxFlow = %f"), p_V->getMaxFlow(0)));      
    maxFlowInitalComputed = true;

auto end = chrono::high_resolution_clock::now();
double time_taken =
    chrono::duration_cast<chrono::nanoseconds>(end - start).count();

time_taken *= 1e-9;
    st2->SetLabel(wxString::Format(wxT("Time MF: %f"), time_taken));
    
    vector <Vertex*> mincut;
    p_V->minCut(mincut);
    st3->SetLabel(wxString::Format(wxT("mincut : %d/%d"), mincut.size(), p_V->getV()));

    int label = 2;
    int position = 0;

    for (unsigned int j = 0; j < m_imageHeight; ++j)
    {
        for (unsigned int k = 0; k < m_imageWidth; ++k)
        {
            position = m_imageWidth * j + k;
            label = p_V->getLabel(position + 1);
            m_myImage[3*position] = label == 0 ? 255 : 0;
            m_myImage[3*position+1] = label == 0 ? 255 : 0;
            m_myImage[3*position+2] = label == 0 ? 255 : 0;
        }
    }

	Refresh(false); // update display
}

//------------------------------------------------------------------------
void MyCanvas::BestSize()
//------------------------------------------------------------------------
{
	SetSize(m_imageWidth, m_imageHeight);	// ideal size for canvas
	GetParent()->SetClientSize(GetSize());	// force the main frame to show the whole canvas
}

//------------------------------------------------------------------------
void MyCanvas::OnPaint(wxPaintEvent& WXUNUSED(event))
//------------------------------------------------------------------------
// update the main window content
{
	wxImage* tempImage;  // the bridge between my image buffer and the bitmap to display

	wxPaintDC dc(this);

	if (m_myImage)
	{
		tempImage = new wxImage(m_imageWidth, m_imageHeight, m_myImage, true); // lend my image buffer...
		m_imageBitmap = wxBitmap(*tempImage, -1); // ...to get the corresponding bitmap
		delete(tempImage);		// buffer not needed any more
		dc.DrawBitmap(m_imageBitmap, 0, 0);
	}
}

//The left mouse button movement event processing
void MyCanvas::OnLeftMove(wxMouseEvent& event)
{
    int verPos = 0;
    if (isLeftDown && event.Dragging()) 
    {
        wxPoint mPos = event.GetPosition();
        wxClientDC dc(this);
        dc.SetPen(*wxWHITE_PEN);
        dc.DrawLine(mPos, penPos);
        dc.SetPen(wxNullPen);
        penPos = mPos;
        st1->SetLabel(wxString::Format(wxT("x: %d"), penPos.x));
        st2->SetLabel(wxString::Format(wxT("y: %d"), penPos.y));
        verPos = penPos.y * m_imageWidth + penPos.x;
        p_V->setLabelVertex(verPos + 1, 0);

        /*ofstream log;
        log.open("obj_test.txt", std::ios_base::app);
        log << penPos.x << " " << penPos.y << " " << int(m_myImage[3 * verPos]) << "\n";
        log.close();*/

        if (!maxFlowInitalComputed) { objDots.push_back(m_myImage[3 * verPos]); }
        else  { p_V->updateEdge(verPos, 0); } // obj 
    }
    else if (isRightDown && event.Dragging()) 
    {
        wxPoint mPos = event.GetPosition();
        wxClientDC dc(this);
        dc.SetPen(*wxBLACK_PEN);
        dc.DrawLine(mPos, penPos);
        dc.SetPen(wxNullPen);
        penPos = mPos;
        st1->SetLabel(wxString::Format(wxT("x: %d"), penPos.x));
        st2->SetLabel(wxString::Format(wxT("y: %d"), penPos.y));
        verPos = penPos.y * m_imageWidth + penPos.x;
        p_V->setLabelVertex(verPos + 1, 1);

        /*ofstream log;
        log.open("bkg_test.txt", std::ios_base::app);
        log << penPos.x << " " << penPos.y << " " << int(m_myImage[3 * verPos]) << "\n";
        log.close();*/

        if (!maxFlowInitalComputed) { bkgDots.push_back(m_myImage[3 * verPos]); }
        else { p_V->updateEdge(verPos, 1); } // bkg 
    }
}


// for objects used white dots
void MyCanvas::OnLeftDown(wxMouseEvent& event)
{
    isLeftDown = true;
    penPos = event.GetPosition();
}

// for background used black dots
void MyCanvas::OnRightDown(wxMouseEvent& event)
{
    isRightDown = true;
    penPos = event.GetPosition();
}

void MyCanvas::OnLeftUp(wxMouseEvent& event)
{
    isLeftDown = false;
    if (maxFlowInitalComputed) ProcessImage(); // for object
}

void MyCanvas::OnRightUp(wxMouseEvent& event)
{
    isRightDown = true;
    if (maxFlowInitalComputed) ProcessImage(); // for object
}




//************************************************************************
//************************************************************************
// Frame class (the main window)
//************************************************************************
//************************************************************************

//------------------------------------------------------------------------
class MyFrame : public wxFrame
	//------------------------------------------------------------------------
{
public:
	MyFrame(const wxString& title, const wxPoint& pos, const wxSize& size);
	MyCanvas* m_canvas; // the canvas inside the main frame
	// Event handlers
protected:
	void OnQuit(wxCommandEvent& event);
	void OnAbout(wxCommandEvent& event);
	void OnOpenImage(wxCommandEvent& WXUNUSED(event));
	void OnSaveImage(wxCommandEvent& WXUNUSED(event));
	void OnProcessImage(wxCommandEvent& WXUNUSED(event));
	void OnClose(wxCloseEvent& event);
	void OnBestSize(wxCommandEvent& WXUNUSED(event));
	bool m_imageLoaded;
	DECLARE_EVENT_TABLE()
};


BEGIN_EVENT_TABLE(MyFrame, wxFrame)
EVT_MENU(ID_LOAD, MyFrame::OnOpenImage)
EVT_MENU(ID_SAVE, MyFrame::OnSaveImage)
EVT_MENU(ID_PROCESS, MyFrame::OnProcessImage)
EVT_MENU(ID_BEST_SIZE, MyFrame::OnBestSize)
EVT_MENU(ID_QUIT, MyFrame::OnQuit)
EVT_MENU(ID_ABOUT, MyFrame::OnAbout)
EVT_CLOSE(MyFrame::OnClose)

END_EVENT_TABLE()

//------------------------------------------------------------------------
MyFrame::MyFrame(const wxString& title, const wxPoint& pos, const wxSize& size)
	: wxFrame((wxFrame*)NULL, -1, title, pos, size)
	//------------------------------------------------------------------------
{
	wxMenu* file_menu = new wxMenu();
	file_menu->Append(ID_LOAD, _T("&Open image..."));
	file_menu->Append(ID_PROCESS, _T("&Process image"));
	file_menu->Append(ID_SAVE, _T("&Save image as..."));
	file_menu->Append(ID_BEST_SIZE, _T("&Best size"));
	file_menu->AppendSeparator();
	file_menu->Append(ID_ABOUT, _T("&About..."));
	file_menu->AppendSeparator();
	file_menu->Append(ID_QUIT, _T("&Exit"));

	wxMenuBar* menuBar = new wxMenuBar();
	menuBar->Append(file_menu, _T("&File"));
	SetMenuBar(menuBar);


	// create the canvas that will manage the image
	m_canvas = new MyCanvas(this, -1, wxDefaultPosition, wxDefaultSize);
	m_imageLoaded = false;
	Centre();
}


//------------------------------------------------------------------------
void MyFrame::OnQuit(wxCommandEvent& WXUNUSED(event))
//------------------------------------------------------------------------
{
	Close(true);
}

//------------------------------------------------------------------------
void MyFrame::OnClose(wxCloseEvent& event)
//------------------------------------------------------------------------
{
	delete m_canvas;
	event.Skip();
}

//------------------------------------------------------------------------
void MyFrame::OnAbout(wxCommandEvent& WXUNUSED(event))
//------------------------------------------------------------------------
{
	wxMessageBox(_T("Image Segmentation \n\n github.com/dies0mbre/"),
		_T(APP_NAME), wxOK | wxICON_INFORMATION);
}

//------------------------------------------------------------------------
void MyFrame::OnProcessImage(wxCommandEvent& WXUNUSED(event))
//------------------------------------------------------------------------
{
	if (m_imageLoaded)
		m_canvas->ProcessImage();
}

//------------------------------------------------------------------------
void MyFrame::OnOpenImage(wxCommandEvent& WXUNUSED(event))
//------------------------------------------------------------------------
{
	wxBitmap bitmap;

	wxString filename = wxFileSelector(_T("Select file"), _T(""), _T(""), _T(""), _T("All files (*.*)|*.*"));
	if (!filename.empty())
	{
		m_canvas->LoadImage(filename);
		m_imageLoaded = true;
	}
}

//------------------------------------------------------------------------
void MyFrame::OnSaveImage(wxCommandEvent& WXUNUSED(event))
//------------------------------------------------------------------------
{
	//	char str[128] = "" ; // proposed file name

	if (!m_imageLoaded)
		return;

	wxString filename = wxFileSelector(_T("Save image as"), _T(""), _T(""), _T("*.bmp"), _T("BMP files (*.bmp)|*.bmp|GIF files (*gif)|*.gif|JPEG files (*jpg)|*.jpg|PNG files (*png)|*.png|TIFF files (*tif)|*.tif|XPM files (*xpm)|*.xpm|All files (*.*)|*.*"), wxFD_SAVE);
	if (!filename.empty())
		m_canvas->SaveImage(filename);
}

//------------------------------------------------------------------------
void MyFrame::OnBestSize(wxCommandEvent& WXUNUSED(event))
//------------------------------------------------------------------------
{
	m_canvas->BestSize();
}


//************************************************************************
//************************************************************************
// Application class
//************************************************************************
//************************************************************************

//------------------------------------------------------------------------
class MyApp : public wxApp
	//------------------------------------------------------------------------
{
	virtual bool OnInit();
};

DECLARE_APP(MyApp)
IMPLEMENT_APP(MyApp) // macro that contains the main() function


//------------------------------------------------------------------------
bool MyApp::OnInit()
//------------------------------------------------------------------------
{
	//support all available image formats
	wxInitAllImageHandlers();

	MyFrame* frame = new MyFrame(_T(APP_NAME), wxDefaultPosition, wxSize(400, 300));
	frame->Show(true);
	SetTopWindow(frame);
	return true;
}
