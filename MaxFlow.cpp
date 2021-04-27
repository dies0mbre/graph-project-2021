// C++ program to implement push-relabel algorithm for
// getting maximum flow of graph
#include <vector>
#include <climits>
#include <iostream>
#include <queue>
using namespace std;

struct Edge
{
    // To store current flow and capacity of edge
    int flow, capacity;
//    int flowResidual,

    // An edge u--->v has start vertex as u and end
    // vertex as v.
    int u, v;

    Edge(int flow, int capacity, int u, int v)
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
    int name, h, e_flow;

    Vertex(int name, int h, int e_flow)
    {
        this->name = name;
        this->h = h;
        this->e_flow = e_flow;
    }
};

// To represent a flow network
class Graph
{
    int V; // No. of vertices
    vector<Vertex> ver;
    vector<Edge> edge;
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
    void updateReverseEdgeFlow(int uAdj, int i, int flow);

public:
    Graph(int V); // Constructor

    // function to add an edge to graph
    void addEdge(int u, int v, int w);
    // Return number in adj[u] for edge u->v
    int findEdge(int u, int v);

    // returns maximum flow from s to t
    int getMaxFlow(int s, int t);

    void PrintCondition();
};

Graph::Graph(int V)
{
    this->V = V;

    // all vertices are initialized with 0 height
    // and 0 excess flow
    for (int i = 0; i < V; i++){
        ver.push_back(Vertex(i, 0, 0));
        adj.push_back(vector<Edge>());
    }
}

int Graph::findEdge(int u, int v)
{
    for (int i=0; i<adj[u].size(); ++i)
    {
        if (adj[u][i].v == v)
        {
            return i;
        }
    }
    return -1;
}

void Graph::addEdge(int u, int v, int capacity)
{
    // flow is initialized with 0 for all edge
//    edge.push_back(Edge(0, capacity, u, v));

    int index = findEdge(u, v);
    if (! (index==-1))
    {
        // means u->v in adj list for u
        cout << "Index!=-1 there is "<< u << "->" <<v<< " in adj list for "<< u << endl;
        adj[u][index].capacity += capacity;
        index = findEdge(v, u);
        adj[v][index].capacity += 0;
    }
    else
    {
        cout << "Index==-1 there is no "<< u << "->" <<v<< " in adj list for "<< u << endl;
        adj[u].push_back(Edge(0, capacity, u, v));
        adj[v].push_back(Edge(0, 0, v, u));
    }
}

void Graph::preprocess(int s)
{
    // Making h of source Vertex equal to no. of vertices
    // Height of other vertices is 0.
    ver[s].h = ver.size();
    int flow;
    Vertex endV(0,0,0);

    //
    for (int i = 0; i < adj[s].size(); i++)
    {
        endV = ver[adj[s][i].v];

        flow = adj[s][i].capacity;
        // Residual capacity is 0
        adj[s][i].flow = 0;
        adj[s][i].capacity = 0; // equal to not having this edge in residual

        // Initialize excess flow for adjacent v
        ver[endV.name].e_flow += flow;
        active.push(&ver[endV.name]);
        cout << "Push in \"active\" node " << ver[endV.name].name << " with excess " << active.back()->e_flow << endl;

        addEdge(endV.name, s, flow);
//        edge.push_back(Edge(0, adj[s][i].flow, edge[i].v, s));
    }
    PrintCondition();
}

// Update reverse flow for flow added on ith Edge for u
void Graph::updateReverseEdgeFlow(int uAdj, int i, int flow)
{
    cout << "Update reverse for " << uAdj << "-->" << adj[uAdj][i].v << "with flow " << flow << endl;
    int u = adj[uAdj][i].v, v = uAdj;

    for (int j = 0; j < adj[u].size(); j++)
    {
        if (adj[u][j].v == v)
        {
            adj[u][j].capacity += flow;
            return;
        }
    }

    // adding reverse Edge in residual graph
    addEdge(u, v, flow);
//    Edge e = Edge(0, flow, u, v);
//    edge.push_back(e);
}

// To push flow from overflowing vertex u
bool Graph::push(int u)
{
    cout << "push for << " << u << endl;
    // Traverse through all edges to find an adjacent (of u)
    // to which flow can be pushed
    for (int i = 0; i < adj[u].size(); i++)
    {
        // if capacity==0 then no flow can be pushed
        if (adj[u][i].capacity == 0)
            continue;

        // Push is only possible if height of adjacent
        // is smaller than height of overflowing vertex
        if (ver[u].h > ver[adj[u][i].v].h && ver[u].e_flow>0)
        {

            // Flow to be pushed is equal to minimum of
            // residual capacity of edge and excess flow.
            int flow = min(adj[u][i].capacity,
                           ver[u].e_flow);
//            cout << "flow " << flow << endl;


            // Reduce excess flow for overflowing vertex
            ver[u].e_flow -= flow;

            // Increase excess flow for adjacent
            ver[adj[u][i].v].e_flow += flow;

//            cout << "push from " << u << "-->" << adj[i][i].v << endl;

            // If v!= source and !=sink
            if (adj[u][i].v && adj[u][i].v!=ver.back().name)
            {
                // and if excess in v was 0, then add in active
                if (!ver[adj[u][i].v].e_flow-flow)
                {
                    active.push(&ver[adj[u][i].v]);
                    cout << "Push in \"active\" node " << adj[u][i].v << " with excess " << active.back()->e_flow << endl;
                }
            }


            // Add residual flow (With capacity 0 and negative
            // flow)
            adj[u][i].capacity -= flow;

            cout << "PUSH " << u <<"-->" << adj[u][i].v << ": " << flow << endl;

            updateReverseEdgeFlow(u, i, flow);
//                PrintCondition();
            return true;
        }
    }
//    PrintCondition();
    return false;
}

// function to relabel vertex u
void Graph::relabel(int u)
{
    cout << "\n RELABEL\n";
    // Initialize minimum height of an adjacent
    int mh = INT_MAX;

    // Pop the vertex from the queue

    cout << "Pop from \"active\" node " << active.front()->name << " with excess " << active.front()->e_flow << endl;
    active.pop();

    // Find the adjacent with minimum height
    for (int i = 0; i < adj[u].size(); i++)
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
    if (! ver[u].e_flow) return;

    active.push(&ver[u]);

    cout << "Push in \"active\" node " << ver[u].name << " with excess " << active.back()->e_flow << endl;
    PrintCondition();
}

// main function for printing maximum flow of graph
int Graph::getMaxFlow(int s, int t)
{
    preprocess(s);

    // loop until none of the Vertex is active
    while (!active.empty())
    {
        int u = active.front()->name;
        if (!push(u))
            relabel(u);
    }

    // ver.back() returns last Vertex, whose
    // e_flow will be final maximum flow
    return ver.back().e_flow;
}


void Graph::PrintCondition()
{
    for (int i = 0; i<ver.size(); ++i)
        cout << "Node " << ver[i].name << " with excess " << ver[i].e_flow << " and distance " << ver[i].h << endl;
    for (int i=0; i<adj.size(); ++i)
    {
        for (int j=0; j<adj[i].size(); ++j)
            cout << "Edge " << adj[i][j].u << "-->" << adj[i][j].v << ", " << adj[i][j].flow << "/"
                 << adj[i][j].capacity << "\n";
    }
    cout << "\n\n";
}


// Driver program to test above functions
int main()
{
    int V = 6;
//    int V = 4;
//    int V = 4;
    Graph g(V);


//    g.addEdge(0, 1, 10000);
//    g.addEdge(0, 2, 10000);
//    g.addEdge(1, 2, 1);
//    g.addEdge(2, 3, 10000);
//    g.addEdge(1, 3, 10000);

//    // Creating above shown flow network
    g.addEdge(0, 1, 16);
    g.addEdge(0, 2, 13);
    g.addEdge(1, 2, 10);
    g.addEdge(2, 1, 4);
    g.addEdge(1, 3, 12);
    g.addEdge(2, 4, 14);
    g.addEdge(3, 2, 9);
    g.addEdge(3, 5, 20);
    g.addEdge(4, 3, 7);
    g.addEdge(4, 5, 4);


//    g.addEdge(0, 1, 10);
//    g.addEdge(0, 2, 10);
//    g.addEdge(1, 3, 8);
//    g.addEdge(2, 3, 5);
//    g.addEdge(1, 2, 5);

    // Initialize source and sink
    int s = 0, t = 5;

// << "Maximum flow is "
    cout << g.getMaxFlow(s, t);
    return 0;
}
