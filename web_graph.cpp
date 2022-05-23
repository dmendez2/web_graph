#include <iostream>
#include<iomanip>
#include <string>
#include<vector>
#include<unordered_map>
#include<algorithm>
using namespace std; 

class graph {
public:
    graph() {
        //Vertex count at initialization is 0
        vertexCount = 0; 
    }

    void InsertEdge(string from, string to) {
        int vertex1; 
        int vertex2; 
        //If we don't yet have the from or the to vertex in our graph then we are adding a vertex and its outdegree to begin with is 0
        //Keep track of the vertex count to use as an "index" 
        //Also push it back into our elements vector to keep track of what we have
        //Initialize our graphMap for this vertex
        if (vertexMap.find(from) == vertexMap.end()) {
            vertexMap[from] = vertexCount; 
            outDegree[vertexCount] = 0; 
            inDegreeMap[vertexCount];
            elements.push_back(from); 
            vertexCount += 1;
        }
        if (vertexMap.find(to) == vertexMap.end()) {
            vertexMap[to] = vertexCount;
            outDegree[vertexCount] = 0;
            inDegreeMap[vertexCount];
            elements.push_back(to); 
            vertexCount += 1;
        }
        //Map to to a vector with a pair that contains the vertex count of to and the name of from (We do to to from to keep track of our indegree)
        //Add an outdegree to our from vertex
        vertex1 = vertexMap[from]; 
        vertex2 = vertexMap[to];
        inDegreeMap[vertex2].push_back(make_pair(vertex1, from)); 
        outDegree[vertex1] += 1; 
    }

    void GetAdjacentList() {
        //Goes through each vertex using our vertex count as an index
        //The only non-zero elements in our adjacency list will be those who have an indegree from another vertex (Each of which have an integer vertex index, and string name)
        //The value at each index will be an integer of vertex count to use as an index and a float of 1/(outdegree of the vertex coming into this vertex)
        //Push this into a vector and then push the vector for this vertex into our adjacency list
        for(int ii = 0; ii < vertexCount; ii += 1){
            vector<pair<int,double>> adjacent;
            for (unsigned int jj = 0; jj < inDegreeMap[ii].size(); jj += 1) {
                double rankval = 1.0f;
                int vertex = inDegreeMap[ii].at(jj).first;
                rankval = rankval / (float)outDegree[vertex];
                adjacent.push_back(make_pair(vertex, rankval));
            }
            adjacencyList.push_back(adjacent); 
        }
    }

    //The initial rank is just 1/(vertex count)
    void GetInitialRankMatrix() {
        for (int ii = 0; ii < vertexCount; ii += 1) {
            double rankval = 1 / (float)vertexCount; 
            rankMatrix.push_back(rankval);
        }
    }
    
    /*A loop that executes the calculation of rank matrix for the total number of power iterations
    We do matrix multiplication between the adjacency list and the rank matrix
    To find the proper columns and rows to multiply we find a value in the adjacency list which has an int index (represents its column) and a float value
    We go to that same index in the rank matrix (represents the row with the same value as the column so for example multiplying column 3 of a row in the adjacency list with row 3 in the rank matrix)
    We now multiply the column and row together and add it to the rank value for this row
    We now continue to multiply columns of the adjacency list at this row with its corresponding row in the rank matrix and continue to add the result to the rank value
    This rank value will be the new rank value for this row of the rank value and pushed back into a temporary matrix
    Continue for the rest of the rows in the adjacency list
    Set the rank matrix equal to the temporary matrix 
    Repeat the process for the rest of the power iterations
    */
    void GetFinalRankMatrix(int power) {
        //int count = 0; 
        for (int ii = 1; ii < power; ii += 1) {
            //PrintRankMatrix(count); 
            //count += 1; 
            vector<double> tempMatrix;
            for (int jj = 0; jj < vertexCount; jj += 1) {
                double rankval = 0; 
                for (unsigned int kk = 0; kk < adjacencyList.at(jj).size(); kk += 1) {
                    int index = adjacencyList.at(jj).at(kk).first; 
                    rankval += adjacencyList.at(jj).at(kk).second * rankMatrix.at(index); 
                }
                tempMatrix.push_back(rankval);
            }
            rankMatrix = tempMatrix;
        }
    }
   
    //Sort the element list lexographically in ascending order so we can print it out corrently
    //Initialize the adjacency list and intitial rank matrix
    //Compute the final rank matrix
    //Print out the rank for each element 
    void PrintRank(int power) {
        sort(elements.begin(), elements.end()); 
        GetAdjacentList(); 
        GetInitialRankMatrix();
        //PrintAdjacencyList();
        GetFinalRankMatrix(power);
       //PrintRankMatrix(power - 1);
        for (unsigned int ii = 0; ii < elements.size(); ii += 1) {
            int index = vertexMap[elements.at(ii)];
            if (ii != elements.size() - 1) {
                cout << fixed;
                cout << elements.at(ii) << " " << setprecision(2) << rankMatrix.at(index) <<  endl;
            }
            else {
                cout << fixed; 
                cout << elements.at(ii) << " " << setprecision(2) << rankMatrix.at(index);
            }
        }
    }

private:
    //Key of a vertex with index int, the value is a vector which stores all the vertices that point to our key vertex (Stores the vertex index int and vertex name string)
    unordered_map<int, vector<pair<int, string>>> inDegreeMap;
    //The key is the vertex index and the value is its outdegree
    unordered_map<int, int> outDegree; 
    //The key is the vertex name and the value is the vertex's index
    unordered_map<string, int> vertexMap; 
    //The adjacency list holds all the non-zero values of the adjacency matrix and the values are calculated by 1/outdegree of the vertex in question at that index
    vector<vector<pair<int,double>>> adjacencyList; 
    //The rank matrix holds the rank for each webpage (Where the index of the rank matrix corresponds to the index of the vertex)
    vector<double> rankMatrix; 
    //Stores all the names of our webpages 
    vector<string> elements; 
    //Counts the number of vertices and used as an index for the vertex
    int vertexCount; 
};

int main()
{
    //Input the number of lines and the number of power iterations
    //Create our graph object
    int numLines; 
    int powerIterations; 
    string from; 
    string to; 
    cin >> numLines; 
    cin >> powerIterations; 
    graph g; 

    //Get the data which consists of a "from" webpage to a "to" webpage
    //Insert an edge for each pair of data points
    for (int ii = 0; ii < numLines; ii += 1) {
        cin >> from; 
        cin >> to; 
        g.InsertEdge(from, to); 
    }
    //Print the Ranks of the webpages
    g.PrintRank(powerIterations); 
}


