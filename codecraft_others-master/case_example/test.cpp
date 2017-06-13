#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <sstream>
#include <set>
using namespace std;

struct EdgeInfo {
    int v, w, c;
};

struct CustomerNodeInfo {
    int u, v, w;
};

const int inf = 0x3f3f3f3f;
const long long infll = 0x3f3f3f3f3f3f3f3f;
long long cost = 0;
set<int> servers;
vector<vector<EdgeInfo> > graph;
vector<CustomerNodeInfo> customer_nodes;
int node_num, edge_num, customer_num, server_cost, need_flow;

bool readtopo(char * filename);
long long judge(char * filename);

int main(int argc, char * argv[]) {
    if (argc != 3) {
        cout << "error: give me the cast file and output file\n";
        return 1;
    }
    if(readtopo(argv[1]))
        cout << judge(argv[2]) << endl;
}

bool readtopo(char * filename) {
    int u, v, c, w;
    fstream fin(filename, ios::in);
    if (!fin.is_open()) {
        cout << "error: fail to open" << filename << endl;
        return false;
    }
    fin >> node_num >> edge_num >> customer_num;
    graph.resize(node_num, vector<EdgeInfo>());
    customer_nodes.resize(customer_num);
    fin >> server_cost;
    for (int i = 0; i < edge_num; ++i) {
        fin >> u >> v >> w >> c;
        graph[u].push_back((EdgeInfo){v, w, c});
        graph[v].push_back((EdgeInfo){u, w, c});
    }
    need_flow = 0;
    for (int i = 0; i < customer_num; ++i) {
        fin >> u >> v >> w;
        customer_nodes[i] = (CustomerNodeInfo){u, v, w};
        need_flow += w;
    }
    fin.close();
    return true;
}

long long judge(char * filename) {
    long long cost = 0;
    int size, n, flow;
    fstream fin(filename, ios::in);
    if (!fin.is_open()) {
        cost = infll;
        return cost;
    }
    vector<int> line;
    string buf;
    getline(fin, buf);
    stringstream ss(buf);
    ss >> size;
    getline(fin, buf);
    for (int i = 0; i < size; ++i) {
        line.clear();
        getline(fin, buf);
        stringstream sss(buf);
        while (sss >> n) {
            line.push_back(n);
        }
        flow = line[line.size() - 1];
        servers.insert(line[0]);
        for (int j = 0; j < line.size() - 3; ++j) {
            bool find = false;
            for (int k = 0; k < graph[line[j]].size(); ++k) {
                if (graph[line[j]][k].v == line[j+1]) {
                    find = true;
                    graph[line[j]][k].w -= flow;
                    cost += graph[line[j]][k].c * flow;
                    if (graph[line[j]][k].w < 0) {
                        cost = infll;
                        fin.close();
                        return cost;
                    }
                    break;
                }
            }
            if (!find) {
                cost = infll;
                fin.close();
                return cost;
            }
        }
        customer_nodes[line[line.size() - 2]].w -= flow;
    }
    for (int i = 0; i < customer_nodes.size(); ++i) {
        if (customer_nodes[i].w > 0) {
            cost = infll;
            fin.close();
            return cost;
        }
    }
    int servers_size = servers.size();
    cost += server_cost * servers_size;
    fin.close();
    return cost;
}
