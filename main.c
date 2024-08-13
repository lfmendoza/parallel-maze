#include <iostream>
#include <vector>
#include <queue>
#include <omp.h>

struct Cell {
    int x, y;
    double cost;
    bool operator<(const Cell& other) const {
        return cost > other.cost;
    }
};

void dijkstra(std::vector<std::vector<int>>& maze, std::vector<std::vector<double>>& dist, int startX, int startY, int endX, int endY) {
    int rows = maze.size();
    int cols = maze[0].size();
    std::priority_queue<Cell> pq;

    pq.push({startX, startY, 0});
    dist[startX][startY] = 0;

    int dx[] = {-1, 1, 0, 0};
    int dy[] = {0, 0, -1, 1};

    while (!pq.empty()) {
        Cell current = pq.top();
        pq.pop();

        #pragma omp parallel for
        for (int i = 0; i < 4; ++i) {
            int newX = current.x + dx[i];
            int newY = current.y + dy[i];

            if (newX >= 0 && newX < rows && newY >= 0 && newY < cols) {
                double newCost = dist[current.x][current.y] + maze[newX][newY];
                if (newCost < dist[newX][newY]) {
                    dist[newX][newY] = newCost;
                    pq.push({newX, newY, newCost});
                }
            }
        }
    }

    std::cout << "Shortest path cost: " << dist[endX][endY] << std::endl;
}

int main() {
    int rows, cols, startX, startY, endX, endY;
    std::cout << "Enter maze dimensions (rows cols): ";
    std::cin >> rows >> cols;

    std::vector<std::vector<int>> maze(rows, std::vector<int>(cols));
    std::vector<std::vector<double>> dist(rows, std::vector<double>(cols, 1e9));

    std::cout << "Enter start coordinates (x y): ";
    std::cin >> startX >> startY;
    std::cout << "Enter end coordinates (x y): ";
    std::cin >> endX >> endY;

    std::cout << "Enter maze matrix:" << std::endl;
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            std::cin >> maze[i][j];
        }
    }

    dijkstra(maze, dist, startX, startY, endX, endY);

    return 0;
}
