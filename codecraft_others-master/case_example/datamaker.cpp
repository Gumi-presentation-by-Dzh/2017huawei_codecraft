#include <iostream>
#include <cstdio>
#include <ctime>
#include <cstdlib>

using namespace std;

const int Cost = 100;
const int maxn = 1000;
const int X = 20; //平均度数 
const int Y = 10; //消费点比例% 
const int vol_maxn = 100;
const int need_maxn = 30;
const int cost_maxn = 10;

struct node {
	int vol, cost;
}; 

struct server {
	int need;
};

node A[maxn][maxn];
server B[maxn];
int main() {
	freopen("data.txt", "w", stdout);
	srand(time(0));
	int node_num = rand() % maxn + 1;
	int server_num = 0;
	int edge_num = 0;
	for (int i = 1; i <= node_num; i++) {
		for (int j = i + 1; j <= node_num; j++) {
			int x = rand() % 10000;
			if (x <= (X * 10000) / node_num) {
				A[i][j].vol = rand() % vol_maxn + 1;
				A[i][j].cost = rand() % cost_maxn + 1;
				++edge_num;
			}
			else A[i][j].vol = 0;
		}
	}
	for (int i = 1; i <= node_num; i++) {
		int y = rand() % 100;
		if (y <= Y) {
			B[i].need = rand() % need_maxn + 1;
			++server_num;
		}
		else B[i].need = 0;
	}
	cout << node_num << ' ' << edge_num << ' ' << server_num << endl;
	cout << endl;
	cout << Cost << endl << endl;
	for (int i = 1; i <= node_num; i++) {
		for (int j = i + 1; j <= node_num; j++)
			if (A[i][j].vol != 0) {
				cout << i-1 << ' ' << j-1 << ' ' << A[i][j].vol << ' ' << A[i][j].cost << endl;
			}
	}
	cout << endl;
	int cnt = 0;
	for (int i = 1; i <= node_num; i++) 
		if (B[i].need != 0){
			cout << cnt++ << ' ' << i-1 << ' ' << B[i].need << endl; 
		}
	cout << endl;
} 
