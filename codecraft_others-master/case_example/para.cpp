#include <stdio.h>
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
#include <fstream>

using namespace std;

const int MAX_V = 1010;
const int MAX_EDGE_NUM = 20000;
const int inf = 0x3f3f3f3f;
const long long infll = 0x3f3f3f3f3f3f3f3f;
pair<long long, int> costs[5][100000];

int kmean_times = 20;
double last_second = 10, pm = 0.1, pc = 0.9, c1 = 2.0, c2 = 2.0, w = 0.8;

struct EdgeInfo {
    int v, w, c;
};

struct CustomerNodeInfo {
    int u, v, w;
};

struct Edge {
    int t, u, c, U;
    Edge *next, *pair;
};

class Particle {
public:
    Particle(int length=0);
    vector<double> v;
    vector<double> v_best;
    vector<double> vp;
    long long cost_best;
    long long cost;
    friend bool operator< (const Particle & p1, const Particle & p2);
};

class LiuXin {
public:
    void readtopo(char * filename);
    void spfa();
    vector<int> kmeans(int k);
    vector<vector<EdgeInfo> > graph;
    vector<CustomerNodeInfo> customer_nodes;
    vector<vector<int> > d;
    bool vis[MAX_V];
    int need_flow, node_num, edge_num, customer_num, server_cost;
};

class JinTao {
public:
    JinTao(LiuXin & lx);
    LiuXin *liuxin;
    void recover();
    void add_edge(int u, int v, int w, int c);
    void add_server(vector<int> &Q);
    int costflow();
    void print_flow(vector<vector<int> > &node, vector<int> &flow);
private:
    int aug(int u, int m);
    bool modlabel();
    Edge epool[MAX_EDGE_NUM * 4 + 10000], *e[MAX_V];
    int psz, s, t;
    int flow, cost, dist, d[MAX_V];
    bool vis[MAX_V];
};

class HGAPSO {
public:
    HGAPSO(JinTao & jt, double pm, double pc, double c1, double c2, double w);
    void initial();
    vector<int> get_best();
    void addone(vector<int> & v);
    void run();
private:
    Particle encode(vector<int> & vini);
    vector<int> decode(vector<double> & v);
    void GA_cross(Particle & s1, Particle & s2);
    void GA_mutation(Particle & s);
    void PSO_update(Particle & s);
    vector<Particle> p;
    Particle gbest;
    int l;
    double GA_pm, GA_pc, PSO_c1, PSO_c2, PSO_w;
    JinTao *jintao;
};


template <class T>
void knuth_shuffle(vector<T> & v);

long long deploy_server(LiuXin & liuxin, JinTao & jintao);

int main(int argc, char * argv[]) {
    char filenames[5][100] = {
        "case0.txt\0",
        "case1.txt\0",
        "case2.txt\0",
        "case3.txt\0",
        "case4.txt\0",
    };
    fstream fout("para.txt", ios::out);
    int _kmean_times = 20;
    double _pm = 0.1, _pc = 0.9, _c1 = 2.0, _c2 = 2.0, _w = 0.8;
    fout << "kmean_times\tpm\tpc\tc1\tc2\tw\tcost\n";
    cout << "kmean_times\tpm\tpc\tc1\tc2\tw\tcost\n";
    long long _cost = infll, cost;
    for (int i = 0; i < 5; ++i) {
        fout << filenames[i] << endl;
        cout << filenames[i] << endl;
        LiuXin liuxin;
        liuxin.readtopo(filenames[i]);
        liuxin.spfa();
        JinTao jintao(liuxin);
        int cnt = 0;
        for (kmean_times = 20; kmean_times <= 21; kmean_times += 5 )
            for (pm = 0.1; pm <= 0.51; pm += 0.1 )
                for (pc = 0.5; pc <= 1.01; pc += 0.1)
                    for (c1 = 1; c1 <= 2.1; c1 += 0.5)
                        for (c2 = 1; c2 <= 2.1; c2 += 0.5)
                            for (w = 0.8; w <= 0.91; w += 0.1) {
                                cost = deploy_server(liuxin, jintao);
                                costs[i][cnt++] = make_pair(cost, 0);
                                cout << kmean_times << "\t" << pm << "\t" << pc << "\t" << c1 << "\t" << c2 << "\t" << w << "\t" << cost << endl;
                            }
        sort(costs[i], costs[i] + cnt);
        int rank = 1;
        costs[i][0].second = 1;
        fout << i << "\t" << 0 << "\t" << 1 << endl;
        for (int j = 1; j < cnt; ++j) {
            if (costs[i][j].first != costs[i][j-1].first)
                ++rank;
            costs[i][j].second = rank;
            fout << i << "\t" << j << "\t" << costs[i][j].second << endl;
        }
    }
    fout.close();
}

long long deploy_server(LiuXin & liuxin, JinTao & jintao) {
    HGAPSO hgapso(jintao, pm, pc, c1, c2, w);
    vector<int> best_server;
    long long best_cost = infll;
    vector<vector<int> > node;
    vector<int> flow;
    int best_index = liuxin.customer_num;
    int block_size = (int)(sqrt(liuxin.customer_num) + 0.1);
    double double_clock_begin = (double)clock() / CLOCKS_PER_SEC, double_clock_end;
    for (int i = 1; i <= liuxin.customer_num; i += block_size) {
        vector<int> server1 = liuxin.kmeans(i), server2;
        best_server = server1;
        jintao.recover();
        jintao.add_server(server1);
        long long best_now = (long long)jintao.costflow() + i * (long long)liuxin.server_cost;
        if (best_now < best_cost) {
            best_cost = best_now;
            best_server = server1;
            best_index = i;
        }
        for (int j = 1; j < kmean_times ; ++j) {
            server2 = liuxin.kmeans(i);
            bool same = true;
            for (int k = 0; k < i; ++k) {
                if (server1[k] != server2[k]) {
                    same = false;
                    break;
                }
            }
            if (same)
                continue;
            jintao.recover();
            jintao.add_server(server2);
            long long cost = (long long)jintao.costflow() + i * (long long)liuxin.server_cost;
            if (cost < best_now) {
                best_now = cost;
                if (best_now < best_cost) {
                    best_cost = best_now;
                    best_server = server2;
                    best_index = i;
                }
            }
            server1 = server2;
        }
    }

    for (int i = max(best_index - block_size, 1); i <= min(best_index + block_size, liuxin.customer_num); i++) {
        vector<int> server1 = liuxin.kmeans(i), server2;
        jintao.recover();
        jintao.add_server(server1);
        hgapso.addone(server1);
        long long best_now = (long long)jintao.costflow() + i * (long long)liuxin.server_cost;
        if (best_now < best_cost) {
            best_cost = best_now;
            best_server = server1;
            best_index = i;
        }
        for (int j = 1; j < kmean_times ; ++j) {
            server2 = liuxin.kmeans(i);
            bool same = true;
            for (int k = 0; k < i; ++k) {
                if (server1[k] != server2[k]) {
                    same = false;
                    break;
                }
            }
            if (same)
                continue;
            jintao.recover();
            jintao.add_server(server2);
            hgapso.addone(server2);
            long long cost = (long long)jintao.costflow() + i * (long long)liuxin.server_cost;
            if (cost < best_now) {
                best_now = cost;
                if (best_now < best_cost) {
                    best_cost = best_now;
                    best_server = server2;
                }
            }
            server1 = server2;
        }
    }
    
    do {
        hgapso.run();
        double_clock_end = (double)clock() / CLOCKS_PER_SEC;
    } while ((double_clock_end -  double_clock_begin) < last_second);

    best_server = hgapso.get_best();
    
    jintao.recover();
    jintao.add_server(best_server);
    best_cost = (long long)jintao.costflow() + best_server.size() * (long long)liuxin.server_cost;
    return best_cost;
}

Particle::Particle(int length): v(length, 0.0), v_best(length, 0.0), vp(length, 0.0), cost_best(infll), cost(infll) {}

bool operator< (const Particle & p1, const Particle & p2) {
    return p1.cost_best == p2.cost_best ? p1.cost < p2.cost : p1.cost_best < p2.cost_best;
}

template <class T>
void knuth_shuffle(vector<T> & v) {
    int i = v.size() - 1, j = 0;
    while (i >= 0) {
        j = rand() % (i + 1);
        swap(v[i], v[j]);
        --i;
    }
}

void LiuXin::readtopo(char * filename) { 
    fstream fin(filename, ios::in);
    if (!fin.is_open()) {
        cout << "error: fail to open" << filename << endl;
        return;
    }
    int u, v, c, w;
    fin >> node_num >> edge_num >> customer_num;
    graph.resize(node_num, vector<EdgeInfo>());
    d.resize(node_num, vector<int>(node_num, inf));
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
}

void LiuXin::spfa() {
    for(int s = 0; s < node_num; ++s) {
        d[s][s] = 0;
        queue<int> Q;
        memset(vis,0,sizeof(vis));
        vis[s] = true;
        Q.push(s);
        while (!Q.empty()) {
            int u = Q.front();
            Q.pop();
            for (unsigned int i = 0; i < graph[u].size(); ++i) {
                int v = graph[u][i].v;
                int dis =  d[s][u] + graph[u][i].c;
                if (dis < d[s][v])
                { 
                    d[s][v] = dis;
                    if (!vis[v])
                    {
                        vis[v] = true;
                        Q.push(v);
                    }
                }
            }
            vis[u] = false;
        }
    }
}

vector<int> LiuXin::kmeans(int k) {
    srand(time(0));
    vector<int> clusters(k);
    int label[customer_num], min_dist, min_index;
    knuth_shuffle(customer_nodes);
    for (int i = 0; i < k; ++i) {
        clusters[i] = customer_nodes[i].v;
    }

    bool update = true;

    while (update) {
        update = false;
        for (int i = 0; i < customer_num; ++i) {
            min_dist= inf;
            min_index = 0;
            for (int j = 0; j < k; ++j) {
                if(d[clusters[j]][customer_nodes[i].v] < min_dist) {
                    min_dist = d[clusters[j]][customer_nodes[i].v];
                    min_index = j;
                }
            }
            if (label[i] != min_index) {
                update = true;
                label[i] = min_index;
            }
        }
        for (int j = 0; j < k; ++j) {
            min_dist = inf;
            min_index = -1;
            for (int l = 0; l < node_num; ++l) {
                int dist = 0;
                for (int i = 0; i < customer_num; ++i) {
                    dist += label[i] == j ? d[l][customer_nodes[i].v] : 0;
                }
                if (dist < min_dist) {
                    min_index = l;
                    min_dist = dist;
                }
            }
            clusters[j] = min_index;
        }
    }
    sort(clusters.begin(), clusters.end());
    return clusters;
}

JinTao::JinTao(LiuXin & lx) {
    liuxin = &lx;
    recover();
}

void JinTao::recover() {
    psz = 0;
    memset(e, 0, sizeof(e));
    for (unsigned int u = 0; u < liuxin->graph.size(); ++u)
        for (unsigned int i = 0; i < liuxin->graph[u].size(); ++i)
            add_edge(u, liuxin->graph[u][i].v, liuxin->graph[u][i].w, liuxin->graph[u][i].c);
    s = liuxin->node_num + liuxin->customer_num;
    t = s + 1;
    for (int i = 0; i < liuxin->customer_num; ++i) {
        add_edge(liuxin->customer_nodes[i].v, liuxin->customer_nodes[i].u + liuxin->node_num, liuxin->customer_nodes[i].w, 0);
        add_edge(liuxin->customer_nodes[i].u + liuxin->node_num, t, liuxin->customer_nodes[i].w, 0);
    }
}

void JinTao::add_edge(int u, int v, int w, int c) {
    Edge *e1 = epool + psz++, *e2 = epool + psz ++;
    *e1 = (Edge){v, w, c, w, e[u], e2}, e[u] = e1;
    *e2 = (Edge){u, 0, -c, 0, e[v], e1}, e[v] = e2;
}


void JinTao::add_server(vector<int> &Q) {
    for (unsigned int i = 0; i < Q.size(); i++) {
        add_edge(s, Q[i], inf, 0);
    }
}

void JinTao::print_flow(vector<vector<int> > & node, vector<int> &flow) {
    while (true) {
        vector<int> Tmp;
        int u = s;
        int S = inf;
        while (u != t) {
            bool flag=false;
            for (Edge *i = e[u]; i; i = i->next) {
                int v = i->t;
                if (i->U > i->u) {
                    S = min(S, i->U - i->u);
                    u = v;
                    flag=true;
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

int JinTao::aug(int u, int m) {
    if (u == t)
        return cost += (long long)dist * m, m;
    int d = m;
    vis[u] = true;
    for (Edge *i = e[u]; i; i = i->next) {
        if (i->u && !i->c && !vis[i->t]) {
            int f = aug(i->t, min(d, i->u));
            i->u -= f;
            i->pair->u += f;
            d -= f;
            if (!d)
                return m;
        }
    }
    return m - d;
}

bool JinTao::modlabel() {
    deque <int > q;
    memset(vis , 0, sizeof(vis));
    memset(d, 0x3f , sizeof(d));
    q.push_back(s);
    d[s] = 0;
    vis[s] = true;
    while (!q.empty ()) {
        int u = q.front ();
        q.pop_front ();
        vis[u] = false;
        for (Edge *i = e[u]; i; i = i->next) {
            int v = i->t;
            if (i->u && d[u] + i->c < d[v]) {
                d[v] = d[u] + i->c;
                if (vis[v])
                    continue;
                vis[v] = true;
                if (q.size () && d[v] < d[q[0]])
                    q.push_front(v);
                else
                    q.push_back(v);
            }
        }
    }
    for (Edge *i = epool; i < epool + psz; ++i) {
        i->c -= d[i->t] - d[i->pair->t];
    }
    dist += d[t];
    return d[t] < inf;
}

int JinTao::costflow() {
    flow = cost = dist = 0;
    while (modlabel ()) {
        int tmpf;
        do {
            memset(vis , 0, sizeof(vis));
            tmpf = aug(s, INT_MAX);
            flow += tmpf;
        } while (tmpf);
    }
    if (flow != liuxin->need_flow)
        cost = inf;
    return cost;
}

HGAPSO::HGAPSO(JinTao & jt, double pm, double pc, double c1, double c2, double w) {
    GA_pm = pm;
    GA_pc = pc;
    PSO_c1 = c1;
    PSO_c2 = c2;
    PSO_w = w;
    jintao = &jt;
    l = jintao->liuxin->node_num;
    gbest = Particle(l);
}

Particle HGAPSO::encode(vector<int> & v) {
    Particle res = Particle(l);
    for (unsigned int i = 0; i < v.size(); ++i) {
        res.v[v[i]] = 1.0;
    }
    return res;
}

vector<int> HGAPSO::decode(vector<double> & v) {
    vector<int> res;
    for (int i = 0; i < l; ++i) {
        if (v[i] > 0.5)
            res.push_back(i);
    }
    return res;
}

void HGAPSO::initial() {
    p.clear();
    gbest = Particle(l);
}

void HGAPSO::addone(vector<int> & v) {
    p.push_back(encode(v));
}

vector<int> HGAPSO::get_best() {
    return decode(gbest.v_best);
}

void HGAPSO::GA_cross(Particle & s1, Particle & s2) {
    int r1 = rand() % l, r2 = rand() % l;
    if (r1 > r2) 
        swap(r1, r2);
    while (r1 < r2) {
        swap(s1.v[r1], s2.v[r1]);
        ++r1;
    }
}

void HGAPSO::GA_mutation(Particle & s) {
    int r1 = rand() % l, r2 = rand() % l;
    if (r1 > r2) 
        swap(r1, r2);
    while (r1 < r2) {
        s.v[r1] += (double)rand() / RAND_MAX - 0.5;
        s.v[r1] = max(0.0, s.v[r1]);
        s.v[r1] = min(s.v[r1], 1.0);
        ++r1;
    }
}

void HGAPSO::PSO_update(Particle & s) {
    for (int i = 0; i < l; ++i) {
        s.vp[i] = PSO_w * s.vp[i] + PSO_c1 * rand() / RAND_MAX * (s.v_best[i] - s.v[i]) + PSO_c2 * rand() / RAND_MAX * (gbest.v_best[i] - s.v[i]); 
        s.v[i] = s.v[i] + s.vp[i];
        s.v[i] = max(0.0, s.v[i]);
        s.v[i] = min(s.v[i], 1.0);
    }
}

void HGAPSO::run() {
    int r1, r2, i = 0, k = p.size(), j = k / 2;
    double r;
    for (int i = 0; i < k; ++i) {
        jintao->recover();
        vector<int> Q = decode(p[i].v);
        jintao->add_server(Q);
        p[i].cost = (long long)jintao->costflow() + Q.size() * (long long)jintao->liuxin->server_cost;
        if (p[i].cost < p[i].cost_best) {
            p[i].v_best = p[i].v;
            p[i].cost_best = p[i].cost;
            if (p[i].cost < gbest.cost_best) {
                gbest.v_best = p[i].v_best;
                gbest.cost_best = p[i].cost_best;
            }
        }
    }
    sort(p.begin(), p.end());
    while (i < j) {
        PSO_update(p[i]);
        ++i;
    }
    while (i < k) {
        r1 = rand() % j;
        r2 = rand() % j;
        p[i] = p[r1] < p[r2] ? p[r1] : p[r2]; 
        ++i;
    }
    i = j;
    while (i < k) {
        r = (double)rand() / RAND_MAX;
        if (r < GA_pc && i + 1 < k)
            GA_cross(p[i], p[i+1]);
        r = (double)rand() / RAND_MAX;
        if (r < GA_pm)
            GA_mutation(p[i]);
        ++i;
    }
    
}
