#include "deploy.h"
#include <stdio.h>
#include <iostream>
#include <deque>
#include <queue>
#include <vector>
#include <limits.h>
#include <iostream>
#include <cstdlib>
#include <algorithm>
#include <random>
#include <cmath>
#include <cstring>
#include <ctime>
#include <utility>
#include <cmath>

using namespace std;

template<typename T> void pvector(vector<T> v);

double last_second = 90 - 1;
double first_second;

void deploy_server(char * topo[MAX_EDGE_NUM], int line_num,char * filename)
{
    char * topo_file;
    Fuck fuck;

    fuck.readtopo(topo, line_num);
    fuck.spfa();
    XJBS xjbs(fuck);

    vector<int> server;
    vector<int> best_server;
    vector<vector<int> > node;
    vector<int> flow;

    int best_index = fuck.customer_num;
    int kmean_times = 2;
    int block_size = (int)(sqrt(fuck.customer_num) + 0.1);

    // @server: current server locations
    // @best_server: the best server locations with lowest cost until now
    fuck.kmeans(1, server);
    // FIX: what work does @addone do?
    xjbs.addone(server);
    // get the initial best server locations and the initial lowest cost
    fuck.kmeans(best_index, best_server);
    long long best_cost = best_index * fuck.server_cost;

    for (int i = best_index - 1; i > 1; i -= block_size) {
        for (int j = 0; j < kmean_times ; ++j) {
            fuck.kmeans(i, server);
            fuck.add_server(server);
            long long cost = fuck.costflow() + i * fuck.server_cost;
            if (cost < best_cost) {
                best_cost = cost;
                best_index = i;
                best_server.swap(server);
            }
        }
    }

    xjbs.addone(best_server);
    block_size = (block_size >> 1);
    int min_index = max(best_index - block_size, 1);
    int max_index = min(best_index + block_size, fuck.customer_num);
    int max_p_size = 10;
    //int max_p_size = 128 / log(best_index << 3);
    //max_p_size &= 0xfffe;
    //++kmean_times;
    for (int i = min_index; i <= max_index; ++i) {
        for (int j = 0; j < kmean_times; ++j) {
            fuck.kmeans(i, server);
            xjbs.addone(server);
        }
    }

    /*  Initialization *
     *  double last_second = 90 - 1;
     *  double first_second;
     */
    last_second -= xjbs.initial(max_p_size);
    first_second = last_second * 0.4;

    while ((double)clock() / CLOCKS_PER_SEC < first_second)
        xjbs.run1();
    xjbs.reproduction();
    while ((double)clock() / CLOCKS_PER_SEC < last_second)
        xjbs.run2();

    xjbs.get_best(best_server);
    fuck.add_server(best_server);
    fuck.costflow();
    fuck.print_flow(node, flow);

    // output the result
    int node_size = node.size();
    topo_file = new char[node_size * MAX_V * 5];
    topo_file[0] = '\0';
    char line[MAX_V * 5];
    char tmp[100];
    sprintf(line, "%d\n\n", node_size);
    strcat(topo_file, line);
    for (int i = 0; i < node_size; ++i) {
        line[0] = '\0';
        int node_size_1 = node[i].size() - 1;
        for (int j = 0; j < node_size_1; ++j) {
            sprintf(tmp, "%d ", node[i][j]);
            strcat(line, tmp);
        }
        sprintf(tmp, "%d ", node[i][node_size_1] - fuck.node_num);
        strcat(line, tmp);
        sprintf(tmp, "%d\n", (int)flow[i]);
        strcat(line, tmp);
        strcat(topo_file, line);
    }
    write_result(topo_file, filename);
    delete []topo_file;
}

Particle::Particle(int length): v(length, 0), v_best(length, 0), vp(length, 0), cost_best(infll), cost(infll) {}

/*
 * @l is fuck's @node_num
*   void XJBS::addone(vector<int> & v) {
*       p.emplace_back(l, v, fuck);
*   }
*/
Particle::Particle(int length, vector<int> & vi, Fuck * & fuck): 
    v(length, 0), v_best(length, 0), vp(length, 0) 
{
    int size = vi.size();
    for (int i = 0; i < size; ++i) {
        v_best[vi[i]] = v[vi[i]] = 1;
    }
    fuck->add_server(vi);
    cost_best = fuck->costflow() + size * fuck->server_cost;
    cost = cost_best;
}

// compare two Particles
bool cmp(const Particle & p1, const Particle & p2) {
    return p1.cost == p2.cost ? p1.cost_best < p2.cost_best : p1.cost < p2.cost;
}

// upset the order of vector
template <class T>
void knuth_shuffle(vector<T> & v) {
    int i = v.size() - 1;
    int j = 0;
    while (i >= 0) {
        j = rand() % (i + 1);
        swap(v[i], v[j]);
        --i;
    }
}

void Fuck::readtopo(char * topo[MAX_EDGE_NUM], int line_num) { 
    int line = 0;
    int u, v, c, w;
    if (line < line_num)
        sscanf(topo[line], "%d %d %d", &node_num, &edge_num, &customer_num);
    s = node_num + customer_num;
    t = s + 1;
    psz = 0;
    memset(e, 0, sizeof(e));
    graph.resize(node_num, vector<EdgeInfo>());
    d.resize(node_num, vector<int>(node_num, inf));

    line += 2;
    sscanf(topo[line], "%lld", &server_cost);
    line += 2;
    
    // read network link infomations 
    for (int i = 0; i < edge_num; ++i, ++line) {
        sscanf(topo[line], "%d %d %d %d", &u, &v, &w, &c);
        graph[u].emplace_back(v, c);
        graph[v].emplace_back(u, c);
        add_edge(u, v, w, c);
        add_edge(v, u, w, c);
    }

    ++line;
    need_flow = 0;

    // read customer node infomations
    for (int i = 0; i < customer_num; ++i, ++line) {
        sscanf(topo[line], "%d %d %d", &u, &v, &w);
        customer_nodes.emplace_back(v, w);
        add_edge(v, u + node_num, w, 0);
        //WARNING: builds edges between consumers and super sink;
        add_edge(u + node_num, t, w, 0);
        need_flow += w;
    }
    psz_tmp = psz;
}

// every pair of vertices's shortest path
// return @d[x][y] represents the total flow cost from x to y
// WARNING: where is @d array's initilization
// @vis[x] maybe means the vertex @x whether in queue @q
void Fuck::spfa() {
    for(int s = 0; s < node_num; ++s) {
        d[s][s] = 0;
        deque<int> q;
        memset(vis,0,sizeof(vis));
        vis[s] = 1;
        q.push_back(s);
        while (!q.empty()) {
            int u = q.front();
            q.pop_front();
            vis[u] = 0;
            for (unsigned int i = 0; i < graph[u].size(); ++i) {
                int v = graph[u][i].v;
                int dis =  d[s][u] + graph[u][i].c;
                if (dis < d[s][v])
                { 
                    d[s][v] = dis;
                    if (!vis[v]) {
                        vis[v] = 1;
                        // WARNING: this condition is unusual
                        // first process the vertex having smaller distance
                        if (q.size () && d[s][v] < d[s][q[0]])
                            q.push_front(v);
                        else
                            q.push_back(v);
                    }
                }
            }
        }
    }
}

// @k: number of communities of k-means
// @clusters: all node id of core in each k-means cluster
void Fuck::kmeans(int k, vector<int> & clusters) {
    clusters.resize(k);
    memset(vis, -1, sizeof(vis));
    vector<vector<int> > kmean_node(k);
    int min_dist;
    int min_index;
    // upset the order of @customer_nodes
    knuth_shuffle(customer_nodes);

    // make cluster's location in the directly network nodes
    for (int i = 0; i < k; ++i) {
        clusters[i] = customer_nodes[i].v;
    }

    while (1) {
        for (int i = 0; i < k; ++i)
            kmean_node[i].clear();
        bool update = 0;

        // first distribute every customer to a proper cluster
        // do one iteration of k-means clustering with all customer_nodes
        // but remember, the direct linked network node represents customer node
        for (int i = 0; i < customer_num; ++i) {
            min_dist = inf;
            min_index = 0;
            for (int j = 0; j < k; ++j) {
                if(d[clusters[j]][customer_nodes[i].v] < min_dist) {
                    min_dist = d[clusters[j]][customer_nodes[i].v];
                    min_index = j;
                }
            }
            if (vis[i] != min_index) {
                update = 1;
                vis[i] = min_index;
            }
            //FIXed: why didn't erase @i from original vector
            //because every iteration @kmean_node will be cleared
            kmean_node[vis[i]].push_back(i);
        }

        if (!update) 
            break;

        for (int j = 0; j < k; ++j) {
            min_dist = inf;
            min_index = 0;
            for (int l = 0; l < node_num; ++l) {
                int dist = 0;

                // calculate the total distance between every network node and
                // each customer node which is clustered into @kmean_node[j]
                for (unsigned int i = 0; i < kmean_node[j].size(); ++i)
                    dist += d[l][ customer_nodes[ kmean_node[j][i] ].v ];

                // if this network node's total distance is smaller than @min_dist,
                // it means this nerwork is a better server location of previous node.
                if (dist < min_dist) {
                    min_index = l;
                    min_dist = dist;
                }
            }
            
            // for each cluster(or partition)(the total number of clusters is @k), 
            // record the network node id which is the best server location.
            clusters[j] = min_index;
        }
    }
}

void Fuck::add_server(vector<int> & v) 
{
    // I don't understand this if statement's usage
    if (psz != psz_tmp) 
    {
        psz = psz_tmp;
        // empty edges stored in @e[s]
        for (Edge *j = e[s]; j; j = j->next) {
            int x = j->t;
            // delete e[x]'s first edge, 
            // because the edge induced to @s is the newest added edge
            e[x] = e[x]->next;
        }
        e[s] = 0;
        Edge *j = epool + psz;
        for (Edge *i = epool; i < j; ++i) {
            i->u = i->U;
            i->c = i->C;
        }
    }

    // I guess id @s is the super source
    // add a edge between super source with node in @v
    // the edge's bandwidth is infinite, and cost is 0
    for (unsigned int i = 0; i < v.size(); ++i) {
        add_edge(s, v[i], inf, 0);
    }
}

// @u and @v are two vertices
// @w is the bandwidth of this edge
// @c is the cost of unit Gbps
// add Edge (u, v) and the reverse edge
void Fuck::add_edge(int u, int v, int w, int c) {
    Edge *e1 = epool + psz++;
    Edge *e2 = epool + psz++;
    // WARNING: edge initilization here
    // t, u, c, U, C, *next, *pair
    *e1 = (Edge){v, w, c, w, c, e[u], e2};
    e[u] = e1;
    *e2 = (Edge){u, 0, -c, 0, -c, e[v], e1};
    e[v] = e2;
}


void Fuck::print_flow(vector<vector<int> > & node, vector<int> &flow) {
    node.clear();
    flow.clear();
    while (1) {
        vector<int> Tmp;
        int u = s;
        int S = inf;
        while (u != t) {
            bool flag=0;
            for (Edge *i = e[u]; i; i = i->next) {
                int v = i->t;
                if (i->U > i->u) {
                    S = min(S, i->U - i->u);
                    u = v;
                    flag = 1;
                    break;
                }
            }
            if (!flag) break;
        }
        if (u != t) break;
        u = s;
        flow.push_back(S);
        while (u != t) {
            for (Edge *i = e[u]; i; i = i->next) {
                int v = i->t;
                if (i->U > i->u) {
                    i->u += S;
                    u = v;
                    break;
                }
            }
            if (u != t) Tmp.push_back(u);
        }
        node.push_back(Tmp);
    }
}

// called in costflow(): tmpf = aug(s, INT_MAX);
// the return value is a volume of flow
// "aug" means "augmenting path"
// this function returns the bandwidth of a augmenting path by start node @u
// it is a recursive algorithm where first start node is @s
int Fuck::aug(int u, int m) {
     //  cost += (long long) dist*m;
     //  return m;
    if (u == t)
        return cost += (long long)dist * m, m;
    int d = m;
    vis[u] = 1;
    for (Edge *i = e[u]; i; i = i->next) {
        // if i->u      != 0 &&  // still have some bandwidth
        //    i->c      == 0 &&  // the cost per Gbps is 0 , 
        //    vis[i->t] == 0     // have not visited it's neighbor
        //FIX: why the cost must be 0?
        if (i->u && !i->c && !vis[i->t]) {
            int f = aug(i->t, min(d, i->u));
            i->u -= f;  // the bandwidth with @f will be used
            i->pair->u += f;
            d -= f;
            if (!d) return m;
        }
    }
    // it means return f1+f2+f3+f4+f5+f6....
    return m - d;
}

// IMPORTANT: @D stores the cost of every node from @s to it
//            @d stores the distances
bool Fuck::modlabel() {

    deque <int> q;
    memset(vis , 0, sizeof(vis));
    memset(D, 0x3f , sizeof(D)); // 0x3f == 63

    q.push_back(s);
    D[s] = 0;
    vis[s] = 1;

    // @vis[x] == 1 means vertex @x is in queue @q
    while (!q.empty ()) {
        int u = q.front ();
        q.pop_front ();
        for (Edge *i = e[u]; i; i = i->next) {
            int v = i->t;
            int dis = D[u] + i->c;
            
            if (i->u && dis < D[v]) {
                D[v] = dis;
                if (!vis[v]) {
                    vis[v] = 1;
                    // add @v into @q
                    if (q.size () && D[v] < D[q[0]]) q.push_front(v);
                    else q.push_back(v);
                }
            }
        }
        vis[u] = 0; // due to q.pop_front();
    }

    for (Edge *i = epool; i < epool + psz; ++i) {
        i->c -= D[i->t] - D[i->pair->t];
    }

    dist += D[t];

    return D[t] < inf; //const int inf = 0x3f3f3f3f; 
}

long long Fuck::costflow() {
    flow = dist = 0;
    cost = 0;
    while (modlabel()) {
        int tmpf;
        do {
            memset(vis , 0, sizeof(vis));
            tmpf = aug(s, INT_MAX);
            flow += tmpf;
        } while (tmpf);
    }
    if (flow != need_flow)
        cost = infll;
    return cost;
}

XJBS::XJBS(Fuck & fk) {
    fuck = &fk;
    l = fuck->node_num;
}

// copy @vd to @vi with a check condition
void XJBS::decode(vector<double> & vd, vector<int> & vi) {
    vi.clear();
    for (int i = 0; i < l; ++i) {
        // FIX: why check condition is 0.5?
        if (vd[i] > 0.5)
            vi.push_back(i);
    }
}

void XJBS::addone(vector<int> & v) {
    p.emplace_back(l, v, fuck);
}

void XJBS::get_best(vector<int> & server) {
    decode(gbest.v_best, server);
}

void XJBS::GA_cross(Particle & s1, Particle & s2) {
    //clock_t t1 = clock();
    int r1 = rand() % l, r2 = rand() % l;
    if (r1 > r2)
        swap(r1, r2);
    while (r1 < r2) {
        swap(s1.v[r1], s2.v[r1]);
        ++r1;
    }

    //cout << "GA_cross:" << (double)(clock()- t1) / CLOCKS_PER_SEC << endl;
}


void XJBS::OBMA(Particle & s) {
    //clock_t t1 = clock();
    int r1, r2;
    do {
        r1 = rand() % l;
    } while (s.v[r1] > 0.5);
    do {
        r2 = rand() % l;
    } while (s.v[r2] < 0.5);
    swap(s.v[r1], s.v[r2]);
    /*
    vector<int> server, unserver;
    for (int i = 0; i < l; ++i) {
        if (s.v[i] > 0.5) {
            server.push_back(i);
        } else {
            unserver.push_back(i);
        }
    }
    int server_size = server.size(), unserver_size = unserver.size(), server_no, unserver_no;
    knuth_shuffle(server);
    knuth_shuffle(unserver);
    if (server_size && unserver_size) {
        swap(s.v[server[0]], s.v[unserver[0]]);
    }*/
    /*
    for (int i = 0, server_index = 0, unserver_index = 0; i < 1 && server_index < server_size && unserver_index < unserver_size; ++i, ++server_index, ++unserver_index) {
        server_no = server[server_index];
        unserver_no = unserver[unserver_index];
        swap(s.v[server_no], s.v[unserver_no]);
        swap(server[server_index], unserver[unserver_index]);
    }*/
    //cout << "OBMA:" << (double)(clock()- t1) / CLOCKS_PER_SEC << endl;
}

inline void XJBS::PSO_update(Particle & s) {
    //clock_t t1 = clock();
    for (int i = 0; i < l; ++i) {
        s.vp[i] = PSO_w * s.vp[i] + PSO_c1 * rand() / RAND_MAX * (s.v_best[i] - s.v[i]) + PSO_c2 * rand() / RAND_MAX * (gbest.v_best[i] - s.v[i]); 
        s.v[i] = (1 / (1 + exp(100*(0.5-(s.v[i] + s.vp[i])))));
    }
    //cout << "PSO:" << (double)(clock()- t1) / CLOCKS_PER_SEC << endl;
}

inline void XJBS::updateone(Particle & s) {
    vector<int> v;
    decode(s.v, v);
    fuck->add_server(v);
    s.cost = fuck->costflow() + v.size() * fuck->server_cost;
    if (s.cost < s.cost_best) {
        s.v_best = s.v;
        s.cost_best = s.cost;
        if (s.cost_best < gbest.cost_best) {
            gbest.v_best = s.v_best;
            gbest.cost_best = s.cost_best;
            cnt = 0;
        }
    }
}

// copy original Particles into @p again
inline void XJBS::reproduction() {
    int k = p.size();
    for (int i = 0; i < k; ++i)
        p.push_back(p[i]);
}

void XJBS::run1() {
    for (int i = 0; i < max_p_size; ++i) {
        OBMA(p[i]);
        updateone(p[i]);
        PSO_update(p[i]);
    }
    if (++cnt > 200) {
        first_second *= 0.5;
        last_second *= 0.5;
        cnt = 0;
    }
    //cout << "run1 " << gbest.cost_best << endl;
}

void XJBS::run2() {
    int i = 0;
    int j = max_p_size >> 1;
    sort(p.begin(), p.end(), cmp);
    for (; i < j; ++i)
        OBMA(p[i]);
    for (; i < max_p_size; ++i)
        PSO_update(p[i]);
    for (i = 0, j = max_p_size - 1; i < j; ++i, --j) {
        GA_cross(p[i], p[j]);
        updateone(p[i]);
        updateone(p[j]);
    }
    if (++cnt > 200) {
        first_second *= 0.5;
        last_second *= 0.5;
        cnt = 0;
    }
    //cout << "run2 " << gbest.cost_best << endl;
}

// size == 10
double XJBS::initial(int size) 
{
    max_p_size = size;
    PSO_c1 = c1;  // const double c1 = 1.0;
    PSO_c2 = c2;  // const double c2 = 1.6;
    PSO_w = w;    // const double w  = 0.9;
    cnt = 0;
    int limit_size = min(max_p_size >> 1, (int)p.size());
    vector<int> v;
    sort(p.begin(), p.end(), cmp);

    gbest = p[0];
    decode(gbest.v_best, v);
    int best_size = v.size() * 0.7;
    if ( best_size == 0 )
        printf("==0\n");
    for (int i = limit_size; i < max_p_size; ++i) {
        fuck->kmeans(best_size, v);
        addone(v);
    }
    clock_t t1 = clock();
    run1();
    clock_t t2 = clock();
    return double(t2 - t1) / CLOCKS_PER_SEC;
}

template<typename T> void
pvector(vector<T> v)
{
    cout << "pvector: ";
    for(int i = 0; i != v.size(); ++i)
        cout << v[i] << " ";
    cout << endl;
}


