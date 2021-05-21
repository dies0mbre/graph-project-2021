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
#include <boost/histogram.hpp>
#include <boost/format.hpp>
using namespace boost::histogram;

using std::vector;
using std::cout;
using std::endl;
using std::queue;
using std::deque;
using std::min;
using std::ifstream;
#define APP_NAME "Image Segmentation"

auto objHist = make_histogram(axis::regular<>(8, 0, 255, "x"));
auto bkgHist = make_histogram(axis::regular<>(8, 0, 255, "x"));

// std::pair <wxPoint, 
vector <int> bkgDots;
vector <int> objDots;


struct Edge
{
    // To store current flow and capacity of edge
    double flow, capacity;

    // An edge u--->v has start vertex as u and end
    // vertex as v.
    int u, v;

    Edge(double flow, double capacity, int u, int v)
    {
        this->flow = flow;
        this->capacity = capacity;
        this->u = u;
        this->v = v;
    }
};

// Represent a Vertex
struct Vertex
{
    int name, h, depth, pixelX, pixelY, label;
    // label = 0-object, 1-bkg, 2-unknown
    double maxB, e_flow;

    Vertex(int name, int h, double e_flow, int x, int y, int label = 2, int depth = -1) :
       name(name),  h(h), e_flow(e_flow), depth(depth), pixelX(x), pixelY(y), label(label) {}
};


void computeHistogram()
{
    objHist.fill(objDots);
    bkgHist.fill(bkgDots);

    /* access to the number of values in the bin    
    std::cout << int(objHist.at(objHist.axis(0).index(90))) << std::endl;
    std::cout << objHist.at(bkgHist.axis(0).index(120)) << std::endl;*/
}


double probabilityValue(int depth, int area, int lambda = 50) // area: 0 - obj, bkg - 1
{
    double probObj = double(objHist.at(objHist.axis(0).index(depth))) / objDots.size();
    double probBkg = double(bkgHist.at(bkgHist.axis(0).index(depth))) / bkgDots.size();
    double sumProb = probObj + probBkg;

    if (sumProb < 1.0e-10 || probObj < 1.0e-10 || probBkg < 1.0e-10) //==0
    {
        return 150.0; // as double max
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


double bValue(Vertex* p, Vertex* q, int sigma = 50, int dist = 1)
{
    return exp(-pow(p->depth - q->depth, 2) / (2 * pow(sigma, 2))) / dist;
}


// To represent a flow network
class Graph
{
    int V; // No. of vertices
    vector<Vertex> ver;
    // adjacency list of residual graph
    vector< vector<Edge> > adj;

    queue <Vertex*> active;

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
    Graph(int V) {
        this->V = V;

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

    int getV()
    {
        return V;
    }

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
            adj[u].push_back(Edge(0, capacity, u, v));
            adj[v].push_back(Edge(0, capacity, v, u));
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
                    ver[verPosition].maxB = *max_element(begin(bValues), end(bValues));
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
                    ver[verPosition].maxB = *max_element(begin(bValues), end(bValues));
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
                    ver[verPosition].maxB = *max_element(begin(bValues), end(bValues));
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
                    ver[verPosition].maxB = *max_element(begin(bValues), end(bValues));
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
                    ver[verPosition].maxB = *max_element(begin(bValues), end(bValues));
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
                    ver[verPosition].maxB = *max_element(begin(bValues), end(bValues));
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
                ver[verPosition].maxB = *max_element(begin(bValues), end(bValues));
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
                ver[verPosition].maxB = *max_element(begin(bValues), end(bValues));
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
                ver[verPosition].maxB = *max_element(begin(bValues), end(bValues));
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
            addEdge(i, source, 1 + ver[i].maxB);
            addEdge(i, sink, 0);
        }
        else if (ver[i].label == 1) // bkg dot
        {
            addEdge(i, source, 0);
            addEdge(i, sink, 1 + ver[i].maxB);
        }
        else // unknown dot
        {
            addEdge(i, source, probabilityValue(ver[i].depth, 1)); // bkg
            addEdge(i, sink, probabilityValue(ver[i].depth, 0));   // obj
        }
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
    int u = adj[uAdj][i].v, v = uAdj;

    for (unsigned int j = 0; j < adj[u].size(); j++)
    {
        if (adj[u][j].v == v)
        {
            adj[u][j].capacity += flow;
            return;
        }
    }
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
    int index, neighbour;
    while (!deq.empty())
    {
        int u = deq.front();
        deq.pop_front();

        for (int i = 0; i < adj[u].size(); ++i)
        {
            neighbour = adj[u][i].v;
            // found residual for u->(v) through v
            index = findEdge(neighbour, u);

            if (index != -1)
            {
                // If the edge presents in residual graph
                if (adj[neighbour][index].capacity != 0 && visited[neighbour] == false)
                {
                    distances[neighbour] = distances[u] + 1;
                    deq.push_back(neighbour);
                    visited[neighbour] = true;
                }
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
    // loop until none of the Vertex is active
    while (!active.empty())
    {
        int u = active.front()->name;
        if (!push(u))
        {
            if (countGB % V) //countGB%(3*V)
            {
                relabel(u);
                countGB += 1;
            }
            else
            {
                //                cout << "else Gb=" << countGB << endl;
                globalRelabeling();
                countGB += 1;
            }
        }
    }

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
	wxPoint penPos;

protected:
	void OnLeftDown(wxMouseEvent& event); // mouse event handler
	void OnRightDown(wxMouseEvent& event);
private:
	int m_imageWidth;
	int m_imageHeight;
	wxBitmap m_imageBitmap;	// used to display the image
	wxImage* m_imageRGB;		// used to load the image
	unsigned char* m_myImage;	// used to process the image

    std::unique_ptr<Graph> p_V; // pointer on representation of image via graph

	void OnPaint(wxPaintEvent& event);
	// void OnLeftMouceClicked(wxMouseEvent& event);

    friend class Graph;

	DECLARE_EVENT_TABLE()
};

BEGIN_EVENT_TABLE(MyCanvas, wxPanel)
EVT_PAINT(MyCanvas::OnPaint)
EVT_LEFT_DOWN(MyCanvas::OnLeftDown)
EVT_RIGHT_DOWN(MyCanvas::OnRightDown)
END_EVENT_TABLE()



//------------------------------------------------------------------------
MyCanvas::MyCanvas(wxWindow* parent, wxWindowID id,
	const wxPoint& pos, const wxSize& size)
	: wxPanel(parent, id, pos, size, wxSUNKEN_BORDER)
	//------------------------------------------------------------------------
{
	st1 = new wxStaticText(this, -1, wxT("X"), wxPoint(10, 10));
	st2 = new wxStaticText(this, -1, wxT("Y"), wxPoint(10, 30));
	m_myImage = NULL;
	m_imageRGB = NULL;
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

	// open image dialog box
	m_imageRGB = new wxImage(fileName, wxBITMAP_TYPE_ANY, -1); // ANY => can load many image formats
	m_imageBitmap = wxBitmap(*m_imageRGB, -1); // ...to get the corresponding bitmap

	m_imageWidth = m_imageRGB->GetWidth();
	m_imageHeight = m_imageRGB->GetHeight();

	m_myImage = (unsigned char*)malloc(m_imageWidth * m_imageHeight * 3);
	memcpy(m_myImage, m_imageRGB->GetData(), m_imageWidth * m_imageHeight * 3);
    st1->SetLabel(wxString::Format(wxT("m_imageWidth = %d"), m_imageWidth));
    st2->SetLabel(wxString::Format(wxT("m_imageHeight = %d"), m_imageHeight));

    /*
    Initialization of graph.
    */
    p_V = std::make_unique<Graph>(m_imageWidth* m_imageHeight + 2); //+2 for sink and source
    p_V->Init(m_imageWidth, m_imageHeight, m_myImage);
    st2->SetLabel(wxString::Format(wxT("getSMTH = %d"), p_V->getV()));

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

    computeHistogram();
    p_V->contactWithTerminals();

    
    st1->SetLabel(wxString::Format(wxT("maxFlow: %f"), p_V->getMaxFlow(0)));

    vector <Vertex*> mincut;
    p_V->minCut(mincut);
    st2->SetLabel(wxString::Format(wxT("mincut : %d/%d"), mincut.size(), p_V->getV()));
    
    /*
    int n, m, a, b, c;
    ifstream infile("C:\\Users\\Asus\\CLionProjects\\study\\graphProject\\MaxFlow-tests\\test_rd07.txt");
    infile >> n >> m;
    Graph g(n);
    while (infile >>a >> b >> c)
    {
        g.addEdge(a - 1, b - 1, c);
    }
    st2->SetLabel(wxString::Format(wxT("maxFlow: %f"), g.getMaxFlow(0)));
    */

	// m_myImage is a monodimentional vector of pixels (RGBRGB...)

    int label = 2;
    int position = 0;

    /*for (unsigned int j = 0; j < i; j += 3)
    {
        label = p_V->ver[].label;
        m_myImage[j] = ver[]
    }*/

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

	/*while (i--)
		m_myImage[i] = 255 - m_myImage[i];*/

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

// for objects used white dots
void MyCanvas::OnLeftDown(wxMouseEvent& event)
{
	if (event.LeftDown()) {
		penPos = event.GetPosition();
		wxClientDC dc(this);
		dc.SetPen(*wxWHITE_PEN);
		dc.DrawPoint(penPos);
		dc.SetPen(wxNullPen);
		st1->SetLabel(wxString::Format(wxT("x: %d"), penPos.x));
		st2->SetLabel(wxString::Format(wxT("y: %d"), penPos.y));

        p_V->setLabelVertex(penPos.y * m_imageWidth + penPos.x + 1, 0);
        objDots.push_back(m_myImage[3 * (penPos.y * m_imageWidth + penPos.x)]);
	}
}

// for background used black dots
void MyCanvas::OnRightDown(wxMouseEvent& event)
{
	if (event.RightDown()) {
		penPos = event.GetPosition();
		wxClientDC dc(this);
		dc.SetPen(*wxBLACK_PEN);
		dc.DrawPoint(penPos);
		dc.SetPen(wxNullPen);
		st1->SetLabel(wxString::Format(wxT("x: %d"), penPos.x));
		st2->SetLabel(wxString::Format(wxT("y: %d"), penPos.y));
        // bkgDots.push_back(std::make_pair(penPos, m_myImage[3 * (penPos.y * m_imageWidth + penPos.x)]));
        p_V->setLabelVertex(penPos.y * m_imageWidth + penPos.x + 1, 1);
		bkgDots.push_back(m_myImage[3 * (penPos.y * m_imageWidth + penPos.x)]);
	}
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

	// Event handlers
protected:
	void OnQuit(wxCommandEvent& event);
	void OnAbout(wxCommandEvent& event);
	void OnOpenImage(wxCommandEvent& WXUNUSED(event));
	void OnSaveImage(wxCommandEvent& WXUNUSED(event));
	void OnProcessImage(wxCommandEvent& WXUNUSED(event));
	void OnClose(wxCloseEvent& event);
	void OnBestSize(wxCommandEvent& WXUNUSED(event));

	MyCanvas* m_canvas; // the canvas inside the main frame
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
	wxMessageBox(_T("How to \n\n- load\n- display\n- process\n- save\n\nan image with wxWidgets (2.8.7)\n\nPascal Bertolino - GIPSA-lab, Grenoble - France\npascal.bertolino@gipsa-lab.fr"),
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
