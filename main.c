#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <time.h>
#include <omp.h>
#include <locale.h>
#include <wchar.h>

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
        {1, 2, 3, 9999, 15, 9999, 13, 14, 15, 9999, 11, 12, 15, 13, 9999, 9},
        {9999, 9999, 5, 9999, 15, 11, 9999, 25, 13, 9999, 10, 9999, 11, 12, 14, 9},
        {25, 22, 15, 9999, 9999, 12, 25, 23, 20, 9999, 12, 9999, 23, 9999, 15, 25},
        {15, 13, 12, 14, 25, 21, 20, 18, 13, 9999, 9999, 9999, 15, 15, 25, 22},
        {15, 12, 11, 9999, 9999, 18, 22, 25, 9999, 9999, 20, 9999, 21, 20, 19, 18},
        {9999, 9999, 9999, 11, 9999, 12, 13, 25, 21, 9999, 9999, 9999, 15, 14, 9999, 25},
        {9999, 9999, 21, 15, 9999, 21, 9999, 11, 12, 9999, 13, 18, 19, 15, 22, 9999},
        {14, 13, 22, 18, 25, 12, 14, 12, 14, 9999, 15, 15, 22, 9999, 15, 9999},
        {12, 15, 11, 9999, 22, 9999, 25, 13, 22, 9999, 9999, 9999, 14, 9999, 12, 14},
        {9999, 9999, 9999, 9999, 12, 15, 11, 9999, 13, 25, 15, 15, 9999, 12, 11, 9999},
        {10, 12, 9999, 13, 22, 15, 15, 9999, 15, 13, 15, 13, 22, 12, 15, 14},
        {22, 25, 15, 15, 15, 9999, 9999, 25, 9999, 13, 15, 25, 13, 9999, 22, 13},
        {11, 9999, 11, 15, 9999, 13, 15, 15, 15, 15, 9999, 25, 13, 14, 12, 11},
        {13, 25, 15, 13, 12, 9999, 9999, 14, 15, 25, 15, 9999, 15, 11, 12, 15},
        {9999, 9999, 9999, 13, 9999, 11, 22, 9999, 15, 9999, 15, 13, 22, 15, 11, 25},
        {13, 12, 11, 9999, 9999, 22, 15, 13, 25, 9999, 11, 13, 12, 13, 22, 9999}
    };

    int maze2[16][16] = {
        {25, 9999, 9999, 14, 11, 9999, 22, 13, 25, 9999, 25, 15, 14, 22, 11, 13},
        {9999, 12, 13, 14, 15, 22, 13, 9999, 15, 13, 14, 9999, 22, 25, 11, 22},
        {11, 14, 13, 22, 11, 15, 15, 9999, 15, 9999, 13, 9999, 15, 12, 9999, 22},
        {25, 25, 22, 15, 15, 15, 11, 12, 9999, 9999, 13, 25, 9999, 11, 25, 22},
        {9999, 9999, 12, 13, 9999, 15, 22, 13, 15, 12, 15, 11, 13, 9999, 25, 13},
        {13, 9999, 15, 13, 15, 9999, 11, 25, 25, 12, 13, 14, 9999, 15, 15, 22},
        {25, 12, 13, 14, 25, 9999, 15, 22, 13, 12, 9999, 11, 25, 9999, 22, 11},
        {13, 25, 9999, 13, 14, 11, 22, 15, 25, 9999, 12, 13, 22, 15, 13, 12},
        {25, 13, 15, 15, 15, 22, 9999, 15, 9999, 9999, 11, 12, 13, 15, 15, 25},
        {15, 9999, 22, 11, 12, 13, 22, 9999, 22, 12, 15, 25, 11, 15, 22, 25},
        {15, 11, 15, 13, 9999, 11, 25, 15, 12, 9999, 13, 11, 9999, 22, 15, 13},
        {12, 13, 25, 13, 12, 15, 22, 13, 9999, 11, 25, 11, 13, 22, 9999, 9999},
        {15, 12, 13, 9999, 15, 9999, 12, 13, 22, 11, 15, 25, 9999, 25, 22, 11},
        {11, 13, 15, 22, 25, 9999, 25, 12, 9999, 9999, 15, 13, 15, 11, 9999, 15},
        {22, 13, 9999, 25, 15, 13, 9999, 15, 12, 15, 13, 9999, 11, 25, 12, 11},
        {15, 9999, 9999, 15, 13, 11, 22, 15, 13, 15, 11, 22, 9999, 15, 15, 12}
    };

    int maze3[16][16] = {
        {11, 15, 15, 25, 12, 9999, 9999, 15, 22, 9999, 12, 15, 25, 15, 11, 9999},
        {9999, 22, 15, 15, 9999, 11, 12, 22, 13, 12, 15, 15, 9999, 15, 12, 15},
        {12, 9999, 9999, 9999, 11, 22, 9999, 12, 15, 9999, 13, 14, 11, 9999, 12, 22},
        {15, 25, 22, 12, 13, 22, 25, 9999, 15, 25, 11, 22, 11, 25, 12, 12},
        {25, 9999, 12, 11, 9999, 9999, 22, 9999, 25, 25, 9999, 11, 15, 15, 13, 14},
        {9999, 25, 12, 22, 9999, 13, 15, 9999, 12, 13, 12, 22, 11, 12, 25, 22},
        {15, 22, 13, 12, 11, 25, 22, 9999, 12, 13, 15, 13, 11, 12, 25, 15},
        {9999, 11, 25, 25, 22, 9999, 12, 25, 12, 15, 13, 11, 12, 9999, 11, 22},
        {25, 22, 9999, 12, 15, 22, 15, 15, 9999, 25, 11, 12, 15, 11, 9999, 9999},
        {11, 15, 15, 13, 22, 9999, 9999, 9999, 22, 25, 11, 12, 9999, 11, 22, 15},
        {15, 11, 25, 22, 25, 9999, 12, 13, 11, 15, 22, 11, 15, 12, 25, 13},
        {12, 13, 15, 15, 25, 9999, 15, 12, 11, 22, 9999, 15, 15, 22, 9999, 15},
        {9999, 22, 9999, 9999, 25, 22, 11, 12, 13, 9999, 12, 13, 15, 25, 12, 9999},
        {22, 11, 25, 12, 15, 22, 15, 15, 11, 9999, 15, 9999, 12, 13, 11, 25},
        {12, 25, 15, 9999, 22, 15, 12, 13, 11, 25, 15, 12, 11, 9999, 25, 22},
        {9999, 12, 15, 12, 9999, 25, 12, 15, 13, 12, 11, 22, 15, 25, 12, 11}
    };

    int maze4[16][16] = {
        {12, 25, 15, 12, 22, 15, 12, 15, 22, 11, 25, 15, 9999, 11, 15, 12},
        {9999, 12, 22, 12, 11, 15, 25, 25, 11, 12, 22, 15, 9999, 15, 11, 25},
        {25, 15, 12, 22, 9999, 25, 15, 25, 11, 9999, 25, 12, 13, 11, 15, 12},
        {11, 22, 12, 11, 15, 12, 15, 11, 22, 25, 25, 11, 12, 22, 11, 15},
        {9999, 25, 12, 15, 9999, 12, 15, 22, 11, 12, 15, 9999, 25, 22, 25, 12},
        {12, 11, 15, 9999, 12, 25, 22, 12, 15, 9999, 12, 15, 11, 25, 22, 11},
        {25, 12, 11, 12, 9999, 15, 11, 25, 12, 15, 25, 22, 15, 9999, 11, 25},
        {11, 12, 25, 22, 15, 15, 12, 9999, 15, 25, 11, 12, 15, 12, 11, 15},
        {15, 9999, 11, 12, 25, 15, 25, 9999, 15, 12, 11, 12, 9999, 15, 25, 12},
        {11, 15, 22, 12, 9999, 25, 12, 11, 22, 25, 15, 11, 9999, 12, 15, 22},
        {25, 12, 15, 11, 15, 12, 25, 12, 15, 9999, 25, 22, 15, 25, 9999, 15},
        {11, 12, 25, 22, 15, 12, 15, 11, 22, 25, 12, 11, 15, 9999, 11, 15},
        {15, 9999, 12, 25, 15, 22, 9999, 25, 12, 9999, 15, 22, 25, 12, 15, 11},
        {22, 12, 25, 11, 15, 22, 12, 11, 25, 12, 15, 25, 22, 11, 12, 9999},
        {12, 15, 9999, 12, 25, 15, 22, 25, 11, 12, 15, 22, 11, 9999, 25, 12},
        {25, 12, 15, 11, 22, 15, 12, 11, 25, 15, 22, 12, 9999, 15, 11, 12}
    };

    int maze5[16][16] = {
        {22, 9999, 15, 12, 11, 15, 9999, 22, 12, 15, 9999, 11, 22, 25, 12, 11},
        {9999, 15, 25, 15, 12, 22, 15, 12, 11, 25, 15, 9999, 22, 11, 15, 12},
        {11, 12, 15, 9999, 12, 11, 25, 15, 11, 22, 15, 25, 11, 9999, 22, 15},
        {25, 22, 11, 15, 25, 22, 12, 9999, 25, 11, 9999, 15, 12, 11, 22, 9999},
        {12, 15, 9999, 12, 11, 22, 15, 12, 25, 11, 9999, 12, 15, 25, 12, 22},
        {9999, 25, 11, 9999, 22, 15, 12, 11, 25, 22, 12, 9999, 22, 11, 25, 15},
        {11, 25, 12, 25, 15, 9999, 25, 15, 11, 22, 12, 11, 9999, 25, 12, 22},
        {12, 15, 9999, 11, 15, 9999, 11, 22, 25, 11, 25, 15, 12, 22, 9999, 25},
        {15, 25, 22, 15, 11, 12, 22, 11, 25, 12, 9999, 15, 25, 11, 12, 9999},
        {22, 12, 9999, 15, 25, 12, 11, 9999, 15, 25, 15, 22, 25, 9999, 11, 15},
        {9999, 12, 15, 11, 12, 22, 15, 9999, 22, 12, 15, 11, 25, 12, 9999, 22},
        {15, 25, 11, 12, 25, 15, 11, 22, 25, 12, 22, 25, 12, 15, 22, 11},
        {12, 11, 22, 15, 12, 25, 15, 11, 9999, 22, 12, 11, 22, 15, 25, 12},
        {22, 12, 15, 9999, 11, 12, 22, 15, 12, 11, 25, 15, 25, 9999, 22, 15},
        {15, 9999, 25, 12, 11, 22, 12, 25, 9999, 22, 15, 11, 12, 25, 11, 15},
        {12, 22, 15, 25, 11, 12, 25, 11, 22, 15, 11, 22, 25, 12, 15, 11}
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
    setlocale(LC_ALL, "");

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
                {
                    clock_t start = clock();
                    dijkstra_sequential(maze, dist, prevX, prevY, startX, startY, endX, endY, N);
                    clock_t end = clock();
                    double time_spent = (double)(end - start) / CLOCKS_PER_SEC;
                    printf("Tiempo de ejecución (Secuencial): %.2f segundos\n", time_spent);
                    print_maze(maze, N, startX, startY, endX, endY);
                }
                break;
            case 2:
                {
                    clock_t start = clock();
                    dijkstra_parallel(maze, dist, prevX, prevY, startX, startY, endX, endY, N);
                    clock_t end = clock();
                    double time_spent = (double)(end - start) / CLOCKS_PER_SEC;
                    printf("Tiempo de ejecución (Paralelo): %.2f segundos\n", time_spent);
                    print_maze(maze, N, startX, startY, endX, endY);
                }
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
