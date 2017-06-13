
#include <vector>
#include <map>
#include <set>
#include <utility>
#include <algorithm>
#include <string>
#include <queue>
#include <deque>
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <ctime>
#include <climits>
#include <unistd.h>
#include "lib_io.h"
#include "lib_time.h"
#include "deploy.h"

using namespace std;

void deploy_server(char * topo[MAX_EDGE_NUM], int line_num,char * filename)
{
    graph g;
    g.init(topo, line_num);
    g.have_superends = 0;
    g.allpairs_sp();

    vector<vid> locations;
    vector<vid> best;
    int best_server_num = 0;

    PSO pso(g);
    g.cluster(1, locations);
    pvector<vid>(locations);
    pso.add_candidate(locations);
    g.cluster(g.c_node_num, best);

    LL best_cost = g.c_node_num * g.server_cost;

    int times = 2;
    costflow cf(g);
    int block_size = (int)(sqrt(g.c_node_num) + 0.1);
    //FIX: the third statement of for loop
    for (int i = g.c_node_num - 1; i > 1; i--)
    {
        for (int j = 0; j < times; ++j)
        {
            g.cluster(i, locations);
            g.make_supersource(locations);
            cf.compute();
            LL cost = cf.getflowcost() + i * g.server_cost;
            if ( cost < best_cost )
            {
                best_cost = cost;
                best_server_num = i;
                best.swap(locations);
            }
        }
    }
    pso.add_candidate(best);
    block_size = block_size >> 1;
    int min_num = max(best_server_num - block_size, 1);
    int max_num = min(best_server_num + block_size, g.c_node_num);
    int max_p_size = 10;

    for(int i = min_num; i <= max_num; ++i)
    {
        for (int j = 0; j < times; ++j)
        {
            g.cluster(i, locations);
            pso.add_candidate(locations);
        }
    }

    pso.last_second = 89 - pso.init(max_p_size);
    pso.first_second = pso.last_second * 0.4;

    while((double)clock() / CLOCKS_PER_SEC < pso.first_second)
        pso.phase(1);
    pso.reproduction();
    while((double)clock() / CLOCKS_PER_SEC < pso.last_second) 
        pso.phase(2);

    pso.get_best(best);
    g.make_supersource(best);
    cf.compute();
    
    vector<vector<vid>> nodes;
    vector<int> flow;
    cf.print_flow(nodes, flow);

    //output
    int node_size = nodes.size();
    
    vid totalnum = g.c_node_num + g.n_node_num;
    //FIX: the @totalnum vs @MAX_V may beyond the bounds of array
    char *topo_file = new char[node_size * (totalnum) * 5];
    topo_file[0] = '\0';
    char line[totalnum * 5];
    char tmp[100];
    sprintf(line, "%d\n\n", node_size);
    strcat(topo_file, line);
    for (int i = 0; i < node_size; ++i)
    {
        line[0] = '\0';
        int node_size_1 = nodes[i].size() - 1;
        for(int j = 0; j < node_size_1; ++j)
        {
            sprintf(tmp, "%d ", nodes[i][j]);
            strcat(line, tmp);
        }
        sprintf(tmp, "%d ", nodes[i][node_size_1] - g.n_node_num);
        strcat(line, tmp);
        sprintf(tmp, "%d\n", (int)flow[i]);
        strcat(line, tmp);
        strcat(topo_file, line);
    }
    write_result(topo_file, filename);
    delete []topo_file;
}

void 
graph::init(char *topo[MAX_EDGE_NUM], int linenum)
{
    int curline = 0;
    sscanf(topo[curline], "%d%d%d", &n_node_num, &edge_num, &c_node_num);
    vid totalnum = n_node_num + c_node_num;

    e.resize(totalnum + 5);
    ssour = totalnum;
    sdest = ssour+1;

    d.resize(n_node_num, vector<int>(n_node_num, INF));
    

    curline += 2;
    sscanf(topo[curline], "%d", &server_cost);
    curline += 2;

    vid sour = 0, dest = 0;
    int tbw = 0, tcost = 0;
    //e.reserve(edge_num + 10);
    for(int i = 0; i < edge_num; ++i)
    {
        sscanf(topo[curline++], "%d%d%d%d", &sour, &dest, &tbw, &tcost);
        addEdge(sour, dest, tbw, tcost, 0);
        addEdge(dest, sour, tbw, tcost, 0);
    }
    
    curline += 1;
    vid cid = 0, nid = 0;
    int tneed = 0;
    consumers.resize(c_node_num);
    for(int i = 0; i < c_node_num; ++i)
    {
        sscanf(topo[curline++], "%d%d%d", &cid, &nid, &tneed);
        //consumers[cid] = consumer(nid, tneed);
        consumers[cid] = consumer(nid, tneed);
        //vid cusid = n_node_num + cid;
        vid cusid = n_node_num + cid;
        addEdge(nid, cusid, tneed, 0, 1);
        addEdge(cusid, sdest, tneed, 0, 1);
        //addEdge(cusid, nid, tneed, 0, 1);
        totalneed += tneed;
    }
}

/* add edge (@_u, @_v), with a pair edge to compute max-flow*/
void
graph::addEdge(vid _u, vid _v, int _bw, int _cost, int flag)
{
    // v, cap, c, re, bw, cost
    e[_u].emplace_back(_v, _bw,  _cost, -1, _bw,  _cost);
    e[_v].emplace_back(_u,   0, -_cost, -1,   0, -_cost);
    e[_u].back().re = e[_v].size()-1;
    e[_v].back().re = e[_u].size()-1;
    if ( flag == 0 )
    {
        e[_u].back().ifre = false;
        e[_v].back().ifre = true;
    }
    if ( flag == 1 )
    {
        e[_u].back().ifre = true;
        e[_v].back().ifre = true;
    }
}

// all-pairs shortest paths just using |V| times spfa algorithm
// only for network nodes whose size is @n_node_num
// fresh the @d with the every edge's cost
void
graph::allpairs_sp()
{
    printf("allpairs_sp()\n");
    vid totalvnum = n_node_num + c_node_num;
    deque<vid> dq;
    vector<bool> if_indq(n_node_num, false);

    //FIXed
    //for(int ver = 0; ver < 2; ++ver)
    for(int ver = 0; ver < n_node_num; ++ver)
    {
        if_indq.clear();
        if_indq.resize(n_node_num, false);
        //printf("ver: %d\n", ver);

        dq.clear();
        dq.push_back(ver);

        if_indq[ver] = true;
        d[ver][ver] = 0;

        //int cnt = 0;
        while(!dq.empty()/* && cnt++ < 3*/)
        {
            vid u = dq.front();
            //printf("u: %d\n", u);
            dq.pop_front();
            
            //for(deque<vid>::iterator it = dq.begin(); it != dq.end(); ++it)
            //    printf("%d ", *it);
            //printf("\n");
            if_indq[u] = false;
            for(int eid = 0; eid < e[u].size(); ++eid)
            {
                if (e[u][eid].ifre == true)
                {
                    //printf("ifre is true\n");
                    continue;
                }
                vid nbr = e[u][eid].v;
                //if (nbr >= n_node_num) continue;
                int newcost = e[u][eid].cost + d[ver][u];
                //printf("u: %d, nbr: %d, ver: %d\n", u, nbr, ver);
                //printf("%d %d\n", e[u][eid].cost, d[ver][nbr]);
                if ( newcost < d[ver][nbr] )
                {
                    d[ver][nbr] = newcost;
                    if ( if_indq[nbr] == false )
                    {
                        if (dq.size() && d[ver][nbr] < d[ver][dq.front()])
                            dq.push_front(nbr);
                        else
                            dq.push_back(nbr);
                        if_indq[nbr] = true;
                    }
                }
            }
        }
    }
    printf("exit allpairs_sp()\n");
}

// do k-means cluster algorithm with many iterations
// @cores: the core of each cluster of k-means
// @k: the required cluster number
void 
graph::cluster(int k, vector<vid> &cores)
{
    printf("enter cluster k:%d\n", k);
    cores.clear();
    cores.resize(k);

    // there are @k communities, every community has some consumers
    // remember: @comm stores the consumer id
    vector<set<vid>> comm(k);   
    set<vid> chosen_id;         // choose a part of network nodes initialize cores
    while(chosen_id.size() < k)
    {
        vid id = rand() % consumers.size();
        chosen_id.insert(id);
    }

    // make cluster's location in the chosen consumers' directly network nodes
    set<vid>::iterator cit = chosen_id.begin();
    for(int i = 0; i < k; ++i)
        cores[i] = consumers[*(cit++)].v;

    vector<int> tags(consumers.size());
    for(int i = 0; i != consumers.size(); ++i)
        tags[i] = (k <= 1) ? 0 : i % (k-1);

    vector<int> mindist(k, INF); // the min totalcost in each cluster
    for(int times = 0; times < MAX_ITERATION; ++times)
    {
        //printf("k-means iteration: %d\n", times);
        /* 1st phase: clustering all consumers into @k clusters */
        bool converaged = true;
        for (int i = 0; i < k; ++i)
            comm[i].clear();
        for (int ci = 0; ci < consumers.size(); ++ci)
        {
            int dist = 0;
            int tmptag = tags[ci];
            for (int i = 0; i < k; ++i){
                if (d[consumers[ci].v][cores[i]] < dist){
                    dist = d[consumers[ci].v][cores[i]];
                    tmptag = i;
                }
            }
            if ( tmptag != tags[ci] )
            {
                converaged = false;
                tags[ci] = tmptag;
            }
            comm[tmptag].insert(ci);
        }
        /* 1st phase end */
        if (converaged)
            return;
        /* 2nd phase: check every network node to find the best core of each cluster */
        for (int i = 0; i < k; ++i)
        {
            for (int ni = 0; ni < n_node_num; ++ni)
            {
                int dist = INF;
                for (set<vid>::iterator it = comm[i].begin(); it != comm[i].end(); ++it)
                {
                    vid co_nid = consumers[*it].v;
                    dist += d[co_nid][ni];
                }
                if (dist < mindist[i]){
                    mindist[i] = dist;
                    cores[i] = ni;
                }
            }
        }
        /* 2nd phase end */
        
    }
}
    
void 
graph::make_supersource(vector<vid> &sours)
{
    for(int i = 0; i < e[ssour].size(); ++i)
    {
        vid u = e[ssour][i].v;
        if( (*(e[u].end()-1)).v != ssour )
            printf("error for not equals @ssour\n");
        e[u].erase(e[u].end()-1);
    }
    e[ssour].clear();
    for (int i = 0; i < sours.size(); ++i)
    {
        addEdge(ssour, sours[i], INF, 0, true);
    }
}

void 
graph::make_supersink()
{
    for (int i = c_node_num + n_node_num; i < ssour; ++i)
    {
        e[i].emplace_back(sdest, INF, 0);
    }
}

costflow::costflow(graph &_g):g(_g), vis(_g.n_node_num + _g.c_node_num + 5, 0)
{
    dist = 0;
}

LL
costflow::getflowcost()
{
    return cost;
}

void
costflow::print_flow(vector<vector<vid>> &nodes, vector<int> &flow)
{
    nodes.clear();
    flow.clear();
    while(true)
    {
        vector<int> tmp;
        int u = g.ssour;
        int S = INF;
        while( u != g.sdest )
        {
            bool flag = 0;
            for ( int i = 0; i < g.e[u].size(); ++i)
            {
                edge &te = g.e[u][i];
                vid v = te.v;
                if (te.bw > te.cap)
                {
                    S = min(S, te.bw - te.cap);
                    u = v;
                    flag = 1;
                    break;
                }
            }
            if ( flag == 0  ) break;
        }
        if ( u != g.sdest ) break;
        u = g.ssour;
        flow.push_back(S);
        while( u != g.sdest )
        {
            for (int i = 0; i < g.e[u].size(); ++i)
            {
                edge &te = g.e[u][i];
                vid v = te.v;
                if ( te.bw > te.cap )
                {
                    te.cap += S;
                    u = v;
                    break;
                }
            }
            if ( u != g.sdest ) tmp.push_back(u);
        }
        nodes.push_back(tmp);
    }
}


/* using zkw algorithm to compute the min-cost max-flow */
LL
costflow::compute()
{
    flow = dist = 0;
    cost = 0;
    int cnt = 0;
    while(modlabel())
    {
        int f = INF;
        while(f){
            resetvis(0);
            f = augmenting(g.ssour, INT_MAX);
            flow += f;
        }
    }
    printf("flow: %d\n", flow);
    if (flow != g.totalneed)
        cost = INF;
    return cost;
}

bool
costflow::modlabel()
{

    D.clear();
    D.resize(g.sdest+2, INF);
    //memset(D, 0x3f, sizeof(D));
    deque<int> dq;
    resetvis(0);

    dq.push_back(g.ssour);
    D[g.ssour] = 0;
    vis[g.ssour] = 1;

    while(!dq.empty()) {
        vid u = dq.front();
        dq.pop_front();
        for ( int i = 0; i < g.e[u].size(); ++i )
        {
            vid v = g.e[u][i].v;
            int dis = D[u] + g.e[u][i].c;
            if (g.e[u][i].cap && dis < D[v])
            {
                D[v] = dis;
                if (vis[v] == 0)
                {
                    vis[v] = 1;
                    if( dq.size() && D[v] < D[dq[0]] ) dq.push_front(v);
                    else dq.push_back(v);
                }
            }
        }
        vis[u] = 0;
    }
    vid totalnum = g.sdest+1;
    for ( int i = 0; i < totalnum; ++i )
    {
        for (int j = 0; j < g.e[i].size(); ++j)
        {
            edge &tmp = g.e[i][j];
            tmp.c -= D[tmp.v] - D[g.e[tmp.v][tmp.re].v];
        }
    }
    dist += D[g.sdest];
    return D[g.sdest] < INF;
}

int
costflow::augmenting(vid v, int m)
{
    if ( v == g.sdest )
    {
        cost += (LL)dist * m;
        return m;
    }
    int d = m;
    vis[v] = 1;
    for (int i = 0; i < g.e[v].size(); ++i)
    {
        edge &te = g.e[v][i];
        if (te.cap && !te.c && !vis[te.v] )
        {
            edge &rete = g.e[te.v][te.re];
            int f = augmenting(te.v, min(d, te.cap));
            te.cap -= f;
            rete.cap += f;
            d -= f;
            if (!d) return m;
        }
    }
    return m - d;
}

void
costflow::resetvis(int value)
{
    vis.clear();
    vis.resize(g.n_node_num+g.c_node_num+5, value);
}

void
costflow::resetD(int value)
{
    //D.clear();
    //D.resize(g.n_node_num + g.c_node_num + 5, value);
}

template<typename T> void
pvector(vector<T> &v, int fd = 1)
{
    string s = "";
    sort(v.begin(), v.end());
    cout << "pvector: ";
    for (typename vector<T>::iterator it = v.begin(); it != v.end(); ++it)
    {
        cout << *it << " ";
    }
    cout << endl;
}

void 
PSO::decode(vector<double> & vd, vector<int> & vi) {
	vi.clear();
    for (int i = 0; i < g.n_node_num; ++i) {
        if (vd[i] > 0.5)
            vi.push_back(i);
    }
}

inline void 
PSO::reproduction() {
    int k = candidates.size();
    for (int i = 0; i < k; ++i)
        candidates.push_back(candidates[i]);
}

void 
PSO::add_candidate(vector<int> & cores) {
    candidates.emplace_back(g.n_node_num, cores, g);		//fuckû�иĶ� 
}

double 
PSO::init(int size) {
	max_size = size;
	c1 = 1.0;
    c2 = 1.6;
    w = 0.9;
    cnt = 0;
	int limit = min(max_size >> 1, (int)candidates.size());	
    vector<int> v;
    sort(candidates.begin(), candidates.end(), cmp);
    gbest = candidates[0];
    
    decode(gbest.v_best, v);
    int best_size = v.size() * 0.7;					//0.7��֪����ô���� 
    
    for (int i = limit; i < max_size; ++i) {
        g.cluster(best_size, v);
        add_candidate(v);
    }
    
    clock_t t1 = clock();
    phase(1);
    clock_t t2 = clock();
    return double(t2 - t1) / CLOCKS_PER_SEC;
}

void 
PSO::GA_cross(Particle & s1, Particle & s2) {
    //clock_t t1 = clock();
    int r1 = rand() % g.n_node_num, r2 = rand() % g.n_node_num;
    if (r1 > r2)
        swap(r1, r2);
    while (r1 < r2) {
        swap(s1.v[r1], s2.v[r1]);
        ++r1;
    }
}

void 
PSO::OBMA(Particle & s) {
    int r1, r2;
    do {
        r1 = rand() % g.n_node_num;
    } while (s.v[r1] > 0.5);
    do {
        r2 = rand() % g.n_node_num;
    } while (s.v[r2] < 0.5);				//0.5��֪��Ϊ�����˶��� 
    swap(s.v[r1], s.v[r2]);
}

inline void 
PSO::PSO_update(Particle & s) {						
    for (int i = 0; i < g.n_node_num; ++i) {
        s.vp[i] = w * s.vp[i] + c1 * rand() / RAND_MAX * (s.v_best[i] - s.v[i]) + c2 * rand() / RAND_MAX * (gbest.v_best[i] - s.v[i]); 
        s.v[i] = (1 / (1 + exp(100*(0.5-(s.v[i] + s.vp[i])))));
    }
}

inline void 
PSO::updateone(Particle & s) {
    vector<int> v;
    decode(s.v, v);
    g.make_supersource(v);
	costflow cf(g); // FIX
    cf.compute();
    s.cost = cf.getflowcost() + v.size() * g.server_cost;
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

void 
PSO::get_best(vector<int> & server) {
    decode(gbest.v_best, server);
}

void 
PSO::phase(int i) {
    switch(i)
    {
    	case 1:{
			for (int i = 0; i < max_size; ++i) {
		        OBMA(candidates[i]);
		        updateone(candidates[i]);
		        PSO_update(candidates[i]);
		    }
		    if (++cnt > 200) {
		        first_second *= 0.5;
		        last_second *= 0.5;					//��֪Ϊ�Σ����ĺ���������ǰ�� 
		        cnt = 0;
		    }
			break;
		}
		case 2:{
			    int i = 0;
			    int j = max_size >> 1;
			    sort(candidates.begin(),candidates.end(), cmp);
			    for (; i < j; ++i)
			        OBMA(candidates[i]);
			    for (; i < max_size; ++i)
			        PSO_update(candidates[i]);
			    for (i = 0, j = max_size - 1; i < j; ++i, --j) {
			        GA_cross(candidates[i], candidates[j]);
			        updateone(candidates[i]);
			        updateone(candidates[j]);
			    }
			    if (++cnt > 200) {
			        first_second *= 0.5;
			        last_second *= 0.5;
			        cnt = 0;
			    }
			break;
		}
	}

}

PSO::PSO(graph & g_in): g(g_in){}

Particle::Particle(int length):
    v(length, 0), v_best(length, 0), vp(length, 0), cost_best(INF), cost(INF)
{}

Particle::Particle(int length, vector<int> & vi, graph& g):
    v(length, 0), v_best(length, 0), vp(length, 0)
{
    int size = vi.size();
    for (int i = 0; i < size; ++i) {
        v_best[vi[i]] = v[vi[i]] = 1;
    }
    costflow cf(g);
    g.make_supersource(vi);
    cf.compute();
    cost_best = cf.getflowcost() + size * g.server_cost;
    cost = cost_best;
}

bool 
cmp(const Particle & p1, const Particle & p2) {
    return p1.cost == p2.cost ? p1.cost_best < p2.cost_best : p1.cost < p2.cost;
}


