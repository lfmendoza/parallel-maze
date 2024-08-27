#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <time.h>
#include <omp.h>

#define WALL_WEIGHT 9999
#define MIN_WEIGHT 1
#define MAX_WEIGHT 25

typedef struct {
    int x, y;
    double cost;
} Cell;

typedef struct {
    Cell* data;
    int size;
    int capacity;
} PriorityQueue;

void initPriorityQueue(PriorityQueue* pq, int capacity) {
    pq->data = (Cell*)malloc(capacity * sizeof(Cell));
    pq->size = 0;
    pq->capacity = capacity;
}

void push(PriorityQueue* pq, Cell cell) {
    if (pq->size == pq->capacity) {
        pq->capacity *= 2;
        pq->data = (Cell*)realloc(pq->data, pq->capacity * sizeof(Cell));
    }
    pq->data[pq->size++] = cell;
    for (int i = pq->size - 1; i > 0 && pq->data[i].cost < pq->data[i-1].cost; --i) {
        Cell temp = pq->data[i];
        pq->data[i] = pq->data[i-1];
        pq->data[i-1] = temp;
    }
}

Cell pop(PriorityQueue* pq) {
    Cell minCell = pq->data[0];
    pq->data[0] = pq->data[--pq->size];
    int i = 0;
    while (i < pq->size) {
        int left = 2 * i + 1;
        int right = 2 * i + 2;
        int smallest = i;

        if (left < pq->size && pq->data[left].cost < pq->data[smallest].cost) {
            smallest = left;
        }
        if (right < pq->size && pq->data[right].cost < pq->data[smallest].cost) {
            smallest = right;
        }
        if (smallest == i) break;

        Cell temp = pq->data[i];
        pq->data[i] = pq->data[smallest];
        pq->data[smallest] = temp;

        i = smallest;
    }
    return minCell;
}

int isEmpty(PriorityQueue* pq) {
    return pq->size == 0;
}

void freePriorityQueue(PriorityQueue* pq) {
    free(pq->data);
}

void print_path(int** prevX, int** prevY, int endX, int endY, int** maze, int startX, int startY) {
    if (prevX[endX][endY] == -1 && prevY[endX][endY] == -1) {
        maze[endX][endY] = '*';
    } else {
        print_path(prevX, prevY, prevX[endX][endY], prevY[endX][endY], maze, startX, startY);
        maze[endX][endY] = '*';
    }
}

// Algoritmo de Dijkstra secuencial
void dijkstra_sequential(int** maze, double** dist, int** prevX, int** prevY, int startX, int startY, int endX, int endY, int N) {
    PriorityQueue pq;
    initPriorityQueue(&pq, 10);

    Cell start = {startX, startY, 0};
    push(&pq, start);
    dist[startX][startY] = 0;

    int dx[] = {-1, 1, 0, 0};
    int dy[] = {0, 0, -1, 1};

    printf("Secuencial - Iniciando recorrido desde (%d, %d)\n", startX, startY);
    while (!isEmpty(&pq)) {
        Cell current = pop(&pq);
        printf("Visitando celda (%d, %d)\n", current.x, current.y);

        for (int i = 0; i < 4; ++i) {
            int newX = current.x + dx[i];
            int newY = current.y + dy[i];

            if (newX >= 0 && newX < N && newY >= 0 && newY < N) {
                double newCost = dist[current.x][current.y] + maze[newX][newY];
                if (newCost < dist[newX][newY]) {
                    dist[newX][newY] = newCost;
                    prevX[newX][newY] = current.x;
                    prevY[newX][newY] = current.y;
                    Cell next = {newX, newY, newCost};
                    push(&pq, next);
                }
            }
        }
    }

    printf("Secuencial - Camino encontrado\n");
    print_path(prevX, prevY, endX, endY, maze, startX, startY);
    printf("Secuencial - Fin del recorrido\n");

    freePriorityQueue(&pq);
}

// Algoritmo de Dijkstra paralelo
void dijkstra_parallel(int** maze, double** dist, int** prevX, int** prevY, int startX, int startY, int endX, int endY, int N) {
    PriorityQueue pq;
    initPriorityQueue(&pq, 10);

    Cell start = {startX, startY, 0};
    push(&pq, start);
    dist[startX][startY] = 0;

    int dx[] = {-1, 1, 0, 0};
    int dy[] = {0, 0, -1, 1};

    printf("Paralelo - Iniciando recorrido desde (%d, %d)\n", startX, startY);
    
    #pragma omp parallel
    {
        while (!isEmpty(&pq)) {
            Cell current;
            
            #pragma omp critical
            {
                if (!isEmpty(&pq)) {
                    current = pop(&pq);
                    printf("Visitando celda (%d, %d)\n", current.x, current.y);
                }
            }
            
            #pragma omp for
            for (int i = 0; i < 4; ++i) {
                int newX = current.x + dx[i];
                int newY = current.y + dy[i];

                if (newX >= 0 && newX < N && newY >= 0 && newY < N) {
                    double newCost = dist[current.x][current.y] + maze[newX][newY];
                    if (newCost < dist[newX][newY]) {
                        dist[newX][newY] = newCost;
                        prevX[newX][newY] = current.x;
                        prevY[newX][newY] = current.y;

                        Cell next = {newX, newY, newCost};
                        
                        #pragma omp critical
                        {
                            push(&pq, next);
                        }
                    }
                }
            }
        }
    }

    printf("Paralelo - Camino encontrado\n");
    print_path(prevX, prevY, endX, endY, maze, startX, startY);
    printf("Paralelo - Fin del recorrido\n");

    freePriorityQueue(&pq);
}

// Función para imprimir el laberinto
void print_maze(int** maze, int N, int startX, int startY, int endX, int endY) {
    maze[startX][startY] = 'A';
    maze[endX][endY] = 'B';

    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            if (maze[i][j] == WALL_WEIGHT) {
                printf("X ");
            } else if (maze[i][j] == 'A' || maze[i][j] == 'B' || maze[i][j] == '*') {
                printf("%c ", maze[i][j]);
            } else {
                printf(". ");
            }
        }
        printf("\n");
    }
}

// Función para elegir el laberinto
void choose_maze(int** maze, int choice) {
    int maze1[16][16] = {
        {1, 2, 1, 4, 25, 9999, 9999, 8, 6, 9, 7, 12, 9999, 9999, 15, 17},
        {9999, 7, 3, 5, 4, 9999, 11, 12, 25, 9999, 10, 13, 15, 9999, 19, 20},
        {5, 4, 9999, 9, 10, 9999, 12, 15, 9999, 9999, 20, 23, 25, 9999, 17, 21},
        {9, 9999, 7, 6, 11, 9999, 19, 23, 9999, 9999, 18, 15, 13, 9999, 9999, 12},
        {15, 14, 9999, 8, 12, 9999, 21, 25, 9999, 20, 17, 11, 9999, 9999, 10, 19},
        {21, 17, 9999, 10, 14, 11, 8, 9999, 22, 21, 25, 15, 7, 9999, 9, 23},
        {25, 9999, 22, 25, 9999, 16, 13, 9999, 9999, 13, 12, 11, 9, 9999, 19, 15},
        {22, 21, 20, 9999, 17, 12, 8, 9999, 25, 20, 18, 9, 9999, 9999, 7, 14},
        {15, 14, 12, 9999, 20, 15, 17, 19, 9999, 18, 25, 16, 9999, 9999, 21, 12},
        {25, 9999, 24, 22, 25, 21, 15, 18, 9999, 17, 20, 22, 11, 9999, 9999, 13},
        {9999, 9999, 9999, 18, 22, 25, 23, 24, 17, 14, 13, 12, 10, 9, 8, 7},
        {19, 20, 21, 22, 23, 9999, 24, 25, 16, 14, 15, 13, 12, 9999, 10, 11},
        {18, 9999, 24, 9999, 9999, 17, 9999, 19, 25, 20, 9999, 14, 13, 12, 11, 10},
        {9999, 9999, 23, 22, 21, 15, 14, 9999, 13, 12, 9999, 9999, 11, 9999, 25, 12},
        {22, 20, 19, 18, 9999, 13, 12, 9999, 9999, 15, 9999, 9999, 16, 17, 18, 9},
        {21, 25, 24, 23, 22, 9999, 18, 19, 20, 9999, 13, 10, 9, 8, 7, 6}
    };

    int maze2[16][16] = {
        {1, 9999, 9999, 2, 5, 9999, 9, 11, 15, 18, 12, 9999, 22, 20, 25, 24},
        {9999, 1, 9999, 6, 7, 9, 13, 15, 16, 19, 9999, 11, 25, 9999, 9999, 23},
        {3, 5, 9999, 8, 10, 15, 9999, 9999, 19, 9999, 9999, 14, 20, 23, 22, 21},
        {7, 9999, 13, 14, 18, 9999, 25, 9999, 22, 9999, 21, 20, 15, 11, 9, 8},
        {9999, 9999, 17, 9999, 9999, 19, 9999, 21, 19, 16, 18, 9999, 15, 12, 10, 9999},
        {9999, 21, 23, 25, 9999, 15, 9999, 20, 14, 17, 9999, 13, 12, 11, 9999, 9999},
        {24, 9999, 9999, 21, 20, 9999, 15, 9999, 22, 9999, 17, 18, 14, 12, 11, 9999},
        {9999, 20, 18, 19, 16, 13, 10, 9999, 21, 22, 9999, 9999, 20, 19, 17, 13},
        {9999, 9999, 9999, 15, 14, 12, 9999, 9999, 9999, 13, 17, 15, 11, 19, 16, 15},
        {20, 17, 13, 16, 14, 12, 9999, 18, 9999, 15, 9999, 9999, 12, 13, 10, 9},
        {9999, 25, 21, 20, 18, 9999, 14, 17, 19, 20, 9999, 16, 11, 10, 7, 9999},
        {24, 23, 9999, 9999, 13, 12, 9999, 9999, 15, 9999, 25, 24, 22, 18, 17, 13},
        {23, 9999, 9999, 9999, 14, 15, 16, 9999, 20, 9999, 9999, 21, 22, 13, 11, 10},
        {22, 9999, 20, 25, 9999, 14, 12, 9999, 19, 9999, 15, 9999, 11, 9, 8, 7},
        {9999, 15, 14, 13, 18, 9999, 9999, 11, 15, 9999, 9999, 14, 13, 12, 16, 10},
        {20, 18, 9999, 22, 25, 23, 21, 9999, 9999, 16, 9999, 9999, 9, 8, 7, 6}
    };

    int maze3[16][16] = {
        {15, 13, 1, 17, 9999, 9999, 15, 20, 19, 25, 13, 17, 21, 9999, 10, 9},
        {16, 14, 10, 13, 11, 14, 9999, 9999, 18, 21, 12, 15, 13, 9999, 9999, 9999},
        {9999, 12, 11, 9999, 9999, 9999, 15, 22, 20, 9999, 25, 23, 9999, 9999, 12, 7},
        {25, 9999, 9999, 17, 18, 21, 25, 22, 9999, 9999, 9999, 9999, 19, 11, 8, 9},
        {15, 13, 22, 21, 20, 18, 15, 14, 9999, 9999, 23, 25, 13, 9999, 10, 9999},
        {25, 9999, 9999, 15, 14, 12, 9999, 9999, 25, 18, 15, 13, 20, 9999, 15, 19},
        {21, 25, 13, 12, 9999, 9999, 14, 13, 20, 15, 18, 9999, 13, 12, 9, 11},
        {15, 14, 9999, 9999, 20, 18, 13, 12, 15, 9999, 13, 14, 9999, 9999, 25, 24},
        {13, 9999, 12, 11, 25, 22, 9999, 20, 21, 13, 14, 17, 25, 9999, 23, 9999},
        {9999, 12, 11, 10, 9, 9999, 9999, 15, 13, 12, 9999, 25, 23, 22, 15, 19},
        {21, 9999, 9999, 9999, 9999, 12, 15, 17, 25, 9999, 13, 12, 14, 9999, 17, 9999},
        {9999, 9999, 9999, 9999, 15, 13, 9999, 14, 25, 9999, 12, 15, 21, 13, 9999, 24},
        {25, 22, 15, 13, 12, 9999, 9999, 9999, 9999, 9999, 20, 17, 14, 13, 12, 11},
        {15, 9999, 9999, 9999, 9999, 25, 22, 13, 17, 25, 13, 9999, 9999, 12, 9999, 24},
        {9999, 9999, 9999, 21, 25, 23, 21, 15, 13, 9999, 14, 17, 12, 9999, 25, 22},
        {21, 9999, 19, 9999, 9999, 9999, 9999, 25, 9999, 9999, 14, 15, 23, 12, 10, 9}
    };

    int maze4[16][16] = {
        {1, 2, 9999, 9999, 25, 9999, 22, 23, 21, 20, 18, 17, 15, 9999, 10, 7},
        {2, 1, 3, 4, 5, 9999, 25, 9999, 19, 18, 9999, 16, 15, 14, 9999, 12},
        {5, 7, 8, 9999, 25, 9999, 18, 15, 16, 17, 9999, 14, 15, 9999, 13, 11},
        {8, 9999, 9999, 9999, 12, 13, 25, 9999, 17, 18, 9999, 9999, 12, 13, 14, 15},
        {7, 6, 5, 4, 9999, 9999, 25, 23, 20, 19, 25, 21, 9999, 17, 16, 15},
        {9999, 12, 11, 15, 13, 25, 20, 18, 9999, 21, 9999, 20, 19, 15, 12, 11},
        {25, 20, 15, 9999, 19, 20, 25, 9999, 14, 13, 9999, 15, 14, 12, 13, 11},
        {24, 9999, 9999, 9999, 15, 9999, 23, 15, 9999, 13, 18, 20, 25, 9999, 12, 9},
        {13, 12, 25, 9999, 24, 9999, 21, 20, 19, 18, 9999, 13, 12, 11, 25, 22},
        {25, 9999, 15, 14, 9999, 15, 25, 9999, 12, 11, 10, 15, 20, 25, 22, 23},
        {12, 9999, 11, 14, 9999, 20, 9999, 19, 9999, 14, 9999, 25, 23, 22, 12, 11},
        {11, 10, 9999, 9999, 15, 17, 9999, 25, 23, 21, 9999, 9999, 19, 17, 15, 14},
        {9999, 9999, 9999, 9999, 15, 14, 9999, 9999, 23, 9999, 9999, 12, 11, 15, 14, 13},
        {20, 19, 18, 17, 9999, 14, 25, 23, 21, 20, 19, 25, 15, 14, 13, 12},
        {12, 11, 10, 9, 9999, 25, 9999, 14, 9999, 15, 25, 9999, 15, 13, 14, 13},
        {9999, 15, 9999, 25, 21, 19, 17, 9999, 15, 25, 23, 25, 19, 14, 13, 12}
    };

    int maze5[16][16] = {
        {25, 23, 9999, 25, 9999, 9999, 23, 25, 9999, 25, 9999, 9999, 22, 9999, 20, 15},
        {24, 9999, 23, 9999, 21, 9999, 20, 9999, 18, 9999, 9999, 15, 18, 9999, 15, 14},
        {23, 22, 9999, 9999, 9999, 9999, 17, 9999, 9999, 25, 23, 22, 25, 9999, 11, 12},
        {22, 21, 20, 9999, 17, 9999, 16, 9999, 25, 23, 22, 25, 15, 9999, 15, 13},
        {9999, 20, 19, 18, 9999, 16, 25, 9999, 22, 15, 14, 13, 12, 15, 17, 12},
        {9999, 9999, 18, 17, 15, 9999, 9999, 13, 12, 9999, 9999, 25, 23, 9999, 22, 10},
        {15, 17, 16, 9999, 15, 14, 13, 9999, 20, 15, 17, 9999, 9999, 9999, 9999, 9},
        {25, 23, 9999, 9999, 15, 14, 9999, 25, 23, 22, 20, 9999, 25, 22, 25, 21},
        {23, 22, 9999, 9999, 15, 17, 9999, 25, 23, 22, 15, 14, 25, 13, 12, 11},
        {25, 9999, 9999, 21, 20, 15, 9999, 25, 23, 9999, 14, 15, 9999, 12, 11, 10},
        {9999, 22, 20, 9999, 18, 17, 9999, 25, 9999, 22, 15, 13, 12, 11, 9999, 15},
        {23, 9999, 19, 18, 15, 16, 25, 23, 15, 17, 19, 20, 21, 9999, 13, 9999},
        {25, 23, 25, 9999, 15, 19, 9999, 21, 9999, 25, 23, 17, 9999, 19, 9999, 15},
        {24, 22, 20, 9999, 23, 25, 15, 25, 9999, 9999, 25, 23, 22, 21, 19, 9999},
        {23, 20, 9999, 15, 9999, 19, 9999, 9999, 25, 9999, 23, 22, 21, 15, 14, 13},
        {22, 21, 25, 9999, 23, 15, 9999, 25, 23, 22, 21, 19, 15, 14, 9999, 9999}
    };

    int* mazes[5] = {&maze1[0][0], &maze2[0][0], &maze3[0][0], &maze4[0][0], &maze5[0][0]};
    int* selected_maze = mazes[choice - 1];

    for (int i = 0; i < 16; ++i) {
        for (int j = 0; j < 16; ++j) {
            maze[i][j] = selected_maze[i * 16 + j];
        }
    }
}

int main() {
    int N = 16, startX = 0, startY = 0, endX = 15, endY = 15;

    int** maze = (int**)malloc(N * sizeof(int*));
    double** dist = (double**)malloc(N * sizeof(double*));
    int** prevX = (int**)malloc(N * sizeof(int*));
    int** prevY = (int**)malloc(N * sizeof(int*));

    for (int i = 0; i < N; ++i) {
        maze[i] = (int*)malloc(N * sizeof(int));
        dist[i] = (double*)malloc(N * sizeof(double));
        prevX[i] = (int*)malloc(N * sizeof(int));
        prevY[i] = (int*)malloc(N * sizeof(int));
    }

    printf("Selecciona el laberinto a recorrer:\n");
    printf("1. Laberinto 1\n");
    printf("2. Laberinto 2\n");
    printf("3. Laberinto 3\n");
    printf("4. Laberinto 4\n");
    printf("5. Laberinto 5\n");
    int option;
    scanf("%d", &option);

    if (option < 1 || option > 5) {
        printf("Opción inválida.\n");
        return 1;
    }

    choose_maze(maze, option);

    int menu_option = 0;
    while (menu_option != 4) {
        printf("\nMenú:\n");
        printf("1. Recorrer laberinto de forma secuencial\n");
        printf("2. Recorrer laberinto de forma paralela\n");
        printf("3. Imprimir laberinto\n");
        printf("4. Salir\n");
        printf("Seleccione una opción: ");
        scanf("%d", &menu_option);

        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                dist[i][j] = DBL_MAX;
                prevX[i][j] = -1;
                prevY[i][j] = -1;
            }
        }

        switch (menu_option) {
            case 1:
                dijkstra_sequential(maze, dist, prevX, prevY, startX, startY, endX, endY, N);
                print_maze(maze, N, startX, startY, endX, endY);
                break;
            case 2:
                dijkstra_parallel(maze, dist, prevX, prevY, startX, startY, endX, endY, N);
                print_maze(maze, N, startX, startY, endX, endY);
                break;
            case 3:
                print_maze(maze, N, startX, startY, endX, endY);
                break;
            case 4:
                printf("Saliendo...\n");
                break;
            default:
                printf("Opción inválida.\n");
                break;
        }
    }

    for (int i = 0; i < N; ++i) {
        free(maze[i]);
        free(dist[i]);
        free(prevX[i]);
        free(prevY[i]);
    }

    free(maze);
    free(dist);
    free(prevX);
    free(prevY);

    return 0;
}
