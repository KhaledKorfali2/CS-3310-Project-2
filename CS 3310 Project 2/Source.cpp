#include<iostream>
#include<fstream>
#include<string>
#include<ctime>
#include<cstdlib>
#include<chrono>
#include<vector>
using namespace std;

class Timer
{
public:
	Timer(const char* name)
		: m_Name(name), m_Stopped(false)
	{
		m_StartTimepoint = std::chrono::high_resolution_clock::now();
	}

	~Timer()
	{
		if (!m_Stopped)
			Stop();
	}

	void Stop()
	{
		auto endTimepoint = std::chrono::high_resolution_clock::now();

		long long start = std::chrono::time_point_cast<std::chrono::microseconds>(m_StartTimepoint).time_since_epoch().count();
		long long end = std::chrono::time_point_cast<std::chrono::microseconds>(endTimepoint).time_since_epoch().count();
		long long duration = end - start;
		//std::cout << m_Name << ": " << duration << "ms\n";


		//output info to a file
		fstream timeTXTFile;
		timeTXTFile.open("Time.txt", ios::app);
		if (timeTXTFile.is_open())
		{
			if (m_Name == "floydWarshall")
			{
				timeTXTFile << m_Name << ": " << duration << std::endl;
				timeTXTFile << "======================================" << std::endl;
				timeTXTFile.close();
			}
			else
			{
				timeTXTFile << m_Name << ": " << duration << std::endl;
				timeTXTFile.close();
			}
		}

		m_Stopped = true;
	}
private:
	const char* m_Name;
	std::chrono::time_point<std::chrono::high_resolution_clock> m_StartTimepoint;
	bool m_Stopped;
};


//**********************************************Function Prototypes**********************************************
//Functions to generate adjacent matricies
double prob();
vector<vector<int>> generateGraph(const int ARRAYSIZE, const double GRAPHDENSITY);

//Functions to display adjacency and cost matricies, as well as shortest path routes
void printSolution(vector<vector<int>>& dist);
void displayMatrix(vector<vector<int>>& matrix);
void printShortestPathsHelper(vector<vector<int>>const& path);

//Functions used for Dijkstra Algorithm
vector<vector<vector<int>>> dijkstra1(vector<vector<int>>& graph);
int minDistance(vector<vector<int>>& dist, vector<vector<bool>>& sptSet, int src);

//Funtions used for Floyed Warshall Algorithm
vector<vector<vector<int>>> floydWarshall2(vector<vector<int>> const& adjMatrix);




//**********************************************constant variables/macros**********************************************
const int ARRAYSIZE = 7;
const double GRAPHDENSITY = 0.35; //0.35 is sparce and 0.75 is dense
const int NUMBEROFITERATIONS = 1;
//const int INF = INT_MAX;
#define INF INT_MAX


int main()
{
	//creates a seed for random numbers
	unsigned seed = time(0);
	srand(seed);
	
	//Creates adjacencyMatrix (i.e. graph) and initializes all values to 0
	vector<vector<int>> adjacencyMatrix(ARRAYSIZE, vector<int>(ARRAYSIZE));
	adjacencyMatrix.assign(ARRAYSIZE, vector<int>(ARRAYSIZE, 0));

	//Creates vector to store both distance matrices and path matrices (distance is in index [0][][] and path is in index [1][][])
	vector<vector<vector<int>>> DijkstraAnswers(2, vector<vector<int>>(ARRAYSIZE, vector<int>(ARRAYSIZE)));
	vector<vector<vector<int>>> FloyedWarshallAnswers(2, vector<vector<int>>(ARRAYSIZE, vector<int>(ARRAYSIZE)));

	
	//Populates Adjacency Matrix with random values
	adjacencyMatrix = generateGraph(ARRAYSIZE, GRAPHDENSITY);



	//test Adjacency Matrix
	/*adjacencyMatrix = 
	{
		{0, 6, INF, 1, INF, INF, INF}, 
		{INF, 0, INF, INF, 1, INF, INF},
		{2, INF, 0, 3, INF, INF, INF},
		{INF, INF, INF, 0, 4, INF, INF},
		{INF, INF, INF, INF, 0, INF, INF},
		{INF, INF, INF, INF, 2, 0, INF},
		{3, INF, 2, INF, INF, 6, 0}
	};*/

	
	////Debug: Stores graph edges in Edges.txt file
	//ofstream outp("Edges.txt");
	//outp << "size:" << ARRAYSIZE << "\n";
	//
	//for (int i = 0; i < ARRAYSIZE; i++)
	//{
	//	for (int j = 0; j < ARRAYSIZE; j++)
	//	{
	//		if (adjacencyMatrix[i][j] != 0)
	//		{
	//			outp << "(" << i << ',' << j << ")" << '\t' << "cost:" << adjacencyMatrix[i][j] << '\t' << '\n';

	//		}
	//	}
	//}


    cout << "Program is running...\n";



	//Runs Dijkstra and FloyedWarshall Algorithms a user specified number of times on the same adjaceny matrix
	//(each iteration uses a random adjacency matrix)
    for (int i = 0; i < NUMBEROFITERATIONS; i++)
    {
		DijkstraAnswers = dijkstra1(adjacencyMatrix);
		FloyedWarshallAnswers = floydWarshall2(adjacencyMatrix);

		//Displays the distance matrices for each algorithm along with the shortest routes for each src (source) vertex
		//Displays Adjacency Matrix
		cout << "----------------------Adjacency Matrix----------------------" << endl;
		displayMatrix(adjacencyMatrix);

		cout << "\n----------------------Dijkstra----------------------" << endl;
		//Prints Adjacency Matrix of the Solution
		printSolution(DijkstraAnswers[0]);
		//Prints all paths taken to reach the shortest distances from each src node
		printShortestPathsHelper(DijkstraAnswers[1]);

		cout << "\n----------------------Floyd Warshall----------------------" << endl;
		//Prints Adjacency Matrix of the Solution
		printSolution(FloyedWarshallAnswers[0]);
		//Prints all paths taken to reach the shortest distances from each src node
		printShortestPathsHelper(FloyedWarshallAnswers[1]);
		
    }

    //take runtimes from Time.txt and average them
    fstream averageTimes;
    fstream runTimes;
    runTimes.open("Time.txt", ios::in);
    averageTimes.open("Averages.txt", ios::app);
    float totalTimeDijkstra = 0;
    float totalTimeFloyedWarshall = 0;
    if (averageTimes.is_open() && runTimes.is_open())
    {
        string line;
        int lineCountKijkstra = 0;
        int lineCountFloyedWarshall = 0;
        while (getline(runTimes, line))
        {
            if (line[0] == 'd')
            {
                totalTimeDijkstra += std::stof(line.erase(0, 9));
                lineCountKijkstra += 1;
            }
            else if (line[0] == 'f')
            {
                totalTimeFloyedWarshall += std::stof(line.erase(0, 14));
                lineCountFloyedWarshall += 1;
            }
        }
        float averageDijkstra = totalTimeDijkstra / lineCountKijkstra;
        averageTimes << "totalTimeDijkstra: " << std::to_string(totalTimeDijkstra) << endl;
        averageTimes << "lineCountKijkstra: " << std::to_string(lineCountKijkstra) << endl;
        averageTimes << "averageDijkstra: " << std::to_string(averageDijkstra) << endl;

        float averageFloyedWarshall = totalTimeFloyedWarshall / lineCountFloyedWarshall;
        averageTimes << "totalTimeFloyedWarshall: " << std::to_string(totalTimeFloyedWarshall) << endl;
        averageTimes << "lineCountFloyedWarshall: " << std::to_string(lineCountFloyedWarshall) << endl;
        averageTimes << "averageFloyedWarshall: " << std::to_string(averageFloyedWarshall) << endl;

		averageTimes << "-----------------------------------------------------------" << endl;

        runTimes.close();
        averageTimes.close();
    }

    cout <<"\nFinished!!!" << endl;


	return 0;
}


//Looks through the verticies that haven't been added to the path tree (i.e. sptSet) and finds the vertex with min distance 
int minDistance(vector<vector<int>>& dist, vector<vector<bool>>& sptSet, int src)
{
	// Initialize min value
	int min = INF;
	int min_index;

	for (int v = 0; v < ARRAYSIZE; v++)
	{
		if (sptSet[src][v] == false && dist[src][v] <= min)
		{
			min = dist[src][v];
			min_index = v;
		}
	}
	return min_index;
}



// Function to print shortest path from source to j using
// prev array
void printPath(vector<vector<int>>& prev, int j, int src)
{
	// Base Case : If j is source
	if (prev[src][j] == -1)
	{
		return;
	}
	printPath(prev, prev[src][j], src);
	cout << j << " ";
}


//Dijkstra Algorithm
vector<vector<vector<int>>> dijkstra1(vector<vector<int>>& graph)
{
	Timer timer("dijkstra");

	vector<vector<vector<int>>> DijkstraAnswers(2, vector<vector<int>>(ARRAYSIZE, vector<int>(ARRAYSIZE)));
	//Stores the shortest distances from src to i
	vector<vector<int>> dist;
	dist.assign(ARRAYSIZE, vector<int>(ARRAYSIZE, INF));

	//Stores which vertextes are included in the shortest path tree
	vector<vector<bool>> sptSet;
	sptSet.assign(ARRAYSIZE, vector<bool>(ARRAYSIZE, false));

	//Stores the shortest path tree
	vector<vector<int>> prev;
	prev.assign(ARRAYSIZE, vector<int>(ARRAYSIZE, -1));

	//Dikstra algorithm is run once for every vertex (i.e. the source vertex changes every iteration)
	for (int src = 0; src < ARRAYSIZE; src++)
	{
		// Distance of source vertex from itself is always 0
		dist[src][src] = 0;

		// Find shortest path for all vertices
		for (int count = 0; count < ARRAYSIZE - 1; count++)
		{
			// Pick the minimum distance vertex from the set of vertices not yet processed. u is always equal to src in first iteration
			int u = minDistance(dist, sptSet, src);

			// Mark the picked vertex as processed
			sptSet[src][u] = true;

			// Update dist value of the adjacent vertices of the picked vertex
			for (int v = 0; v < ARRAYSIZE; v++)
				//Relaxes edges
				if (!sptSet[src][v] && (graph[u][v] != INF) && (dist[src][u] != INF) && (dist[src][u] + graph[u][v] < dist[src][v]))
				{
					prev[src][v] = u;
					dist[src][v] = dist[src][u] + graph[u][v];
				}
		}
	}

	////Debug: Prints Adjacency Matrix of the Solution
	//printSolution(dist);
	////Debug: Prints all paths taken to reach the shortest distances from each src node
	//printShortestPathsHelper(prev);

	DijkstraAnswers[0] = dist;
	DijkstraAnswers[1] = prev;
	

	////Debug: Print prev matrix
	//displayMatrix(prev);
	////Debug: Prints Adjacency Matrix of the Solution
	//printSolution(DijkstraAnswers[0]);
	////Debug: Prints all paths taken to reach the shortest distances from each src node
	//printShortestPathsHelper(DijkstraAnswers[1]);
	return DijkstraAnswers;
}


// A utility function to print solution 
void printSolution(vector<vector<int>>& dist)
{
	cout << "The following matrix shows the shortest distances  between every pair of vertices" << endl;
	for (int i = 0; i < ARRAYSIZE; i++)
	{
		for (int j = 0; j < ARRAYSIZE; j++)
		{
			if (dist[i][j] == INF)
			{
				cout << "INF" << "	 ";
			}
			else
			{
				cout << dist[i][j] << "	 ";
			}
		}
		cout << endl;
	}
}


// Recursive function to print path of given vertex u from source vertex v
void printPath(vector<vector<int>>const& path, int v, int u)
{
	if (path[v][u] == v) {
		return;
	}
	printPath(path, v, path[v][u]);
	cout << path[v][u] << " --> ";
}

// Function to print the shortest cost path information between
// all pairs of vertices
void printShortestPathsHelper(vector<vector<int>>const& path)
{
	for (int v = 0; v < ARRAYSIZE; v++)
	{
		for (int u = 0; u < ARRAYSIZE; u++)
		{
			if (u != v && path[v][u] != -1)
			{
				cout << "The shortest path from " << v << " to " << u << ": " << v << " --> ";
				printPath(path, v, u);
				cout << u << endl;
			}
		}
	}
}

// Function to run the Floyd–Warshall algorithm
vector<vector<vector<int>>> floydWarshall2(vector<vector<int>> const& adjMatrix)
{
	Timer timer("floydWarshall");

	vector<vector<vector<int>>> FloyedWarshallAnswers(2, vector<vector<int>>(ARRAYSIZE, vector<int>(ARRAYSIZE)));

	//Stores the shortest distance from vertex v to u 
	vector<vector<int>> cost(ARRAYSIZE, vector<int>(ARRAYSIZE));
	vector<vector<int>> prev(ARRAYSIZE, vector<int>(ARRAYSIZE));

	// initialize cost[] and path[]
	for (int v = 0; v < ARRAYSIZE; v++)
	{
		for (int u = 0; u < ARRAYSIZE; u++)
		{
			// initially, cost would be the same as the weight of the edge
			cost[v][u] = adjMatrix[v][u];

			if (v == u) 
			{
				prev[v][u] = 0;
			}
			else if (cost[v][u] != INF)
			{
				prev[v][u] = v;
			}
			else
			{
				prev[v][u] = -1;
			}
		}
	}

	// run Floyd–Warshall
	for (int k = 0; k < ARRAYSIZE; k++)
	{
		for (int v = 0; v < ARRAYSIZE; v++)
		{
			for (int u = 0; u < ARRAYSIZE; u++)
			{
				//Updates cost and path if taking vertex k results in a shorter distance being traveresed from v to u
				if (cost[v][k] != INF && cost[k][u] != INF && cost[v][k] + cost[k][u] < cost[v][u])
				{
					cost[v][u] = cost[v][k] + cost[k][u];
					prev[v][u] = prev[k][u];
				}
			}
		}
	}

	////Debug: Prints the shortest distance matrix
	//printSolution(cost);

	////Debug: Print the shortest path between all pairs of vertices
	//printShortestPathsHelper(path);

	////Debug: Print prev matrix
	//displayMatrix(prev);

	FloyedWarshallAnswers[0] = cost;
	FloyedWarshallAnswers[1] = prev;

	////Debug: Prints the shortest distance matrix
	//printSolution(FloyedWarshallAnswers[0]);

	////Debug: Print the shortest path between all pairs of vertices
	//printShortestPathsHelper(FloyedWarshallAnswers[1]);

	return FloyedWarshallAnswers;
}


//Generates a random number between 0 and 1
double prob()
{
	return (static_cast<double>(rand()) / RAND_MAX);
}

//Displays matrices
void displayMatrix(vector<vector<int>>& matrix)
{
	for (int i = 0; i < ARRAYSIZE; i++)
	{
		for (int j = 0; j < ARRAYSIZE; j++)
		{
			if (matrix[i][j] == INF)
			{
				cout << "INF\t";
			}
			else
			{
				cout << matrix[i][j] << "\t";
			}
		}
		cout << "\n";
	}
}

//Generates random adjacency matrices
vector<vector<int>> generateGraph(const int ARRAYSIZE, const double GRAPHDENSITY)
{
	//Stores which edges should be assined a cost
	vector<vector<bool>> graph;
	//Stores the cost associated with each edge
	vector<vector<int>> cost;

	//Initailzes graph and cost vectors to a value of 0
	graph.assign(ARRAYSIZE, vector<bool>(ARRAYSIZE, 0));
	cost.assign(ARRAYSIZE, vector<int>(ARRAYSIZE, 0));

	for (int i = 0; i < ARRAYSIZE; ++i)
	{
		for (int j = 0; j < ARRAYSIZE; ++j)
		{
			if (i == j) // No vertex can have an edge from itself to itself
			{
				graph[i][j] = false;
				//cout << graph[i][j];
			}
			else // Randomly decides which vertecies should get edges based on the graph density (higher density = more edges, and vice versa)
			{
				graph[i][j] = (prob() < GRAPHDENSITY);
				//cout << graph[i][j];
			}
		}
	}

	for (int i = 0; i < ARRAYSIZE; ++i)
	{
		for (int j = 0; j < ARRAYSIZE; ++j)
		{
			if (i != j) // No vertex can have an edge from itself to itself
			{
				//Randomally assings cost values to the indicies which were previously designated as true
				if (graph[i][j] == true) 
				{
					cost[i][j] = prob() * 30;
				}
				else //There are no edges between these vertices
				{
					cost[i][j] = INF;
				}
			}
			else // Makes the distance to travel from a vertex to itself 0
			{
				cost[i][j] = 0;
			}
		}
	}

	return cost;
}
