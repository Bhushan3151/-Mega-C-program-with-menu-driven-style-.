/* All_In_One.c
   Mega demo:
   - Global SIGINT handler (Ctrl+C) -> timestamp + exit
   - Polling demo (q to quit)
   - Number system: even/odd, prime, factorial, Fibonacci, GCD, LCM
   - Algorithms: linear/binary search, many sorting algorithms
   - IPC demo: pipe, SysV message queue, shared memory with semaphore
   - Device-driver simulation: ring buffer, writer thread, reader (blocking) + SIGUSR1 handler
   - Bitwise operations
   - Data Structures: Queue, Singly LL, Doubly LL, Circular LL, Heap, Binary Tree, Graph
*/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/ipc.h>
#include <sys/msg.h>
#include <sys/shm.h>
#include <pthread.h>
#include <signal.h>
#include <time.h>
#include <termios.h>
#include <fcntl.h>
#include <semaphore.h>
#include <errno.h>
#include <sys/wait.h>



/*======================= Color_Code ==========================*/
#define BLINK "\033[5m"
#define RESET   "\033[0m"
#define BLACK   "\033[0;30m"
#define RED     "\033[0;31m"
#define GREEN   "\033[0;32m"
#define YELLOW  "\033[0;33m"
#define BLUE    "\033[0;34m"
#define PURPLE  "\033[0;35m"
#define CYAN    "\033[0;36m"
#define WHITE   "\033[0;37m"

#define BOLD_BLACK   "\033[1;30m"
#define BOLD_RED     "\033[1;31m"
#define BOLD_GREEN   "\033[1;32m"
#define BOLD_YELLOW  "\033[1;33m"
#define BOLD_BLUE    "\033[1;34m"
#define BOLD_PURPLE  "\033[1;35m"
#define BOLD_CYAN    "\033[1;36m"
#define BOLD_WHITE   "\033[1;37m"

#define UNDERLINE_BLACK   "\033[4;30m"
#define UNDERLINE_RED     "\033[4;31m"
#define UNDERLINE_GREEN   "\033[4;32m"
#define UNDERLINE_YELLOW  "\033[4;33m"
#define UNDERLINE_BLUE    "\033[4;34m"
#define UNDERLINE_PURPLE  "\033[4;35m"
#define UNDERLINE_CYAN    "\033[4;36m"
#define UNDERLINE_WHITE   "\033[4;37m"

/* -------------------- Global Interrupt (Ctrl+C) -------------------- */
void interrupt_handler(int signum) {
    time_t now;
    time(&now);
    char *timestamp = ctime(&now);

    printf("\n\n*** Interrupt received! Signal number: %d ***\n", signum);
    printf("Interrupt occurred at: %s", timestamp);
    // perform any cleanup here if needed
    exit(0);  // Exit after interrupt
}

/* -------------------- kbhit (polling) -------------------- */
int kbhit() {
    struct termios oldt, newt;
    int ch;
    int oldf;

    // Save old terminal settings
    tcgetattr(STDIN_FILENO, &oldt);
    newt = oldt;

    // Disable canonical mode and echo
    newt.c_lflag &= ~(ICANON | ECHO);
    tcsetattr(STDIN_FILENO, TCSANOW, &newt);

    // Set non-blocking input
    oldf = fcntl(STDIN_FILENO, F_GETFL, 0);
    fcntl(STDIN_FILENO, F_SETFL, oldf | O_NONBLOCK);

    ch = getchar();

    // Restore settings
    tcsetattr(STDIN_FILENO, TCSANOW, &oldt);
    fcntl(STDIN_FILENO, F_SETFL, oldf);

    if (ch != EOF) {
        ungetc(ch, stdin);
        return 1;
    }

    return 0;
}

/* -------------------- Utility functions -------------------- */
void inputArray(int arr[], int n){
    printf("Enter %d elements (space/newline separated):\n", n);
    for(int i = 0; i < n; i++){
        scanf("%d", &arr[i]);
    }
}

void printArray(int arr[], int n){
    for (int i = 0; i < n; i++){
        printf("%d ",arr[i]);
    }
    printf("\n");
}

void swap(int *x, int *y){
    int temp = *x;
    *x = *y;
    *y = temp;
}

/* -------------------- Sorting algorithms -------------------- */
void insertionSort(int arr[], int n){
    for(int i = 1; i < n ; ++i){
        int key = arr[i];
        int j = i - 1;

        while(j >= 0 && arr[j] > key){
            arr[j + 1] = arr[j];
            j = j - 1;
        }
        arr[j + 1] = key;
    }
}

void selectionSort(int arr[], int n){
    for(int i = 0; i < n - 1; i++){
        int min_idx = i;
        for (int j = i + 1; j < n; j++){
            if(arr[j] < arr[min_idx]){
                min_idx = j;
            }
        }
        swap(&arr[i], &arr[min_idx]);
    }
}

void bubbleSort(int arr[], int n){
    int i, j;
    bool swapped;
    for (i = 0; i < n - 1; i++) {
        swapped = false;
        for (j = 0; j < n - i - 1; j++) {
            if (arr[j] > arr[j + 1]) {
                swap(&arr[j], &arr[j + 1]);
                swapped = true;
            }
        }
        if (!swapped) break;
    }
}

int partition(int arr[], int low, int high) {
    int pivot = arr[high];
    int i = low - 1;
    for (int j = low; j <= high - 1; j++) {
        if (arr[j] < pivot) {
            i++;
            swap(&arr[i], &arr[j]);
        }
    }
    swap(&arr[i + 1], &arr[high]);
    return i + 1;
}

void quickSort(int arr[], int low, int high) {
    if (low < high) {
        int pi = partition(arr, low, high);
        quickSort(arr, low, pi - 1);
        quickSort(arr, pi + 1, high);
    }
}

void merge(int arr[], int l, int m, int r) {
    int i, j, k;
    int n1 = m - l + 1;
    int n2 = r - m;
    int *L = malloc(n1 * sizeof(int));
    int *R = malloc(n2 * sizeof(int));
    for (i = 0; i < n1; i++) L[i] = arr[l + i];
    for (j = 0; j < n2; j++) R[j] = arr[m + 1 + j];
    i = 0; j = 0; k = l;
    while (i < n1 && j < n2) {
        if (L[i] <= R[j]) arr[k++] = L[i++];
        else arr[k++] = R[j++];
    }
    while (i < n1) arr[k++] = L[i++];
    while (j < n2) arr[k++] = R[j++];
    free(L); free(R);
}

void mergeSort(int arr[], int l, int r) {
    if (l < r) {
        int m = l + (r - l) / 2;
        mergeSort(arr, l, m);
        mergeSort(arr, m + 1, r);
        merge(arr, l, m, r);
    }
}

void heapify(int arr[], int n, int i) {
    int largest = i;
    int l = 2 * i + 1;
    int r = 2 * i + 2;
    if (l < n && arr[l] > arr[largest]) largest = l;
    if (r < n && arr[r] > arr[largest]) largest = r;
    if (largest != i) {
        swap(&arr[i], &arr[largest]);
        heapify(arr, n, largest);
    }
}

void heapSort(int arr[], int n) {
    for (int i = n / 2 - 1; i >= 0; i--) heapify(arr, n, i);
    for (int i = n - 1; i > 0; i--) {
        swap(&arr[0], &arr[i]);
        heapify(arr, i, 0);
    }
}

void countSort(int inputArray[], int N) {
    int M = 0;
    for (int i = 0; i < N; i++) if (inputArray[i] > M) M = inputArray[i];
    int* countArray = (int*)calloc(M + 1, sizeof(int));
    for (int i = 0; i < N; i++) countArray[inputArray[i]]++;
    for (int i = 1; i <= M; i++) countArray[i] += countArray[i - 1];
    int* outputArray = (int*)malloc(N * sizeof(int));
    for (int i = N - 1; i >= 0; i--) {
        outputArray[countArray[inputArray[i]] - 1] = inputArray[i];
        countArray[inputArray[i]]--;
    }
    for (int i = 0; i < N; i++) inputArray[i] = outputArray[i];
    free(countArray);
    free(outputArray);
}

void bucketSort(float arr[], int n) {
    int i, j, k;
    int bucket_count = 10;
    float buckets[bucket_count][n]; // max size n per bucket
    int bucket_sizes[bucket_count];
    for (i = 0; i < bucket_count; i++) bucket_sizes[i] = 0;

    for (i = 0; i < n; i++) {
        int idx = bucket_count * arr[i]; // assuming 0 <= arr[i] < 1
        if (idx < 0) idx = 0;
        if (idx >= bucket_count) idx = bucket_count - 1;
        buckets[idx][bucket_sizes[idx]++] = arr[i];
    }
    for (i = 0; i < bucket_count; i++)
        for (j = 0; j < bucket_sizes[i]-1; j++)
            for (k = 0; k < bucket_sizes[i]-1-j; k++)
                if (buckets[i][k] > buckets[i][k+1]) {
                    float temp = buckets[i][k];
                    buckets[i][k] = buckets[i][k+1];
                    buckets[i][k+1] = temp;
                }

    int idx = 0;
    for (i = 0; i < bucket_count; i++)
        for (j = 0; j < bucket_sizes[i]; j++)
            arr[idx++] = buckets[i][j];
}

void cocktailSort(int arr[], int n) {
    int swapped = 1;
    int start = 0, end = n-1;
    while (swapped) {
        swapped = 0;
        for (int i = start; i < end; ++i)
            if (arr[i] > arr[i+1]) { swap(&arr[i], &arr[i+1]); swapped = 1; }
        if (!swapped) break;
        swapped = 0;
        end--;
        for (int i = end-1; i >= start; --i)
            if (arr[i] > arr[i+1]) { swap(&arr[i], &arr[i+1]); swapped = 1; }
        start++;
    }
}

int getMax(int arr[], int n) {
    int mx = arr[0];
    for (int i = 1; i < n; i++)
        if (arr[i] > mx) mx = arr[i];
    return mx;
}

void countSortRadix(int arr[], int n, int exp) {
    int *output = malloc(n * sizeof(int));
    int count[10] = {0};

    for (int i = 0; i < n; i++)
        count[(arr[i]/exp)%10]++;
    for (int i = 1; i < 10; i++)
        count[i] += count[i-1];
    for (int i = n-1; i >=0; i--) {
        output[count[(arr[i]/exp)%10]-1] = arr[i];
        count[(arr[i]/exp)%10]--;
    }
    for (int i = 0; i < n; i++) arr[i] = output[i];
    free(output);
}

void radixSort(int arr[], int n) {
    int m = getMax(arr, n);
    for (int exp = 1; m/exp > 0; exp *= 10)
        countSortRadix(arr, n, exp);
}

void shellSort(int arr[], int n) {
    for (int gap = n/2; gap > 0; gap /= 2) {
        for (int i = gap; i < n; i++) {
            int temp = arr[i];
            int j;
            for (j = i; j >= gap && arr[j - gap] > temp; j -= gap)
                arr[j] = arr[j - gap];
            arr[j] = temp;
        }
    }
}

/* -------------------- Searching algorithms -------------------- */
int binary_search(int arr[], int left, int right, int key) {
    while (left <= right) {
        int mid = left + (right - left) / 2;
        if (arr[mid] == key) return mid;
        if (arr[mid] < key) left = mid + 1;
        else right = mid - 1;
    }
    return -1;
}

int linear_search(int arr1[], int N, int X) {
    for (int i = 0; i < N; i++)
        if (arr1[i] == X) return i;
    return -1;
}

int recursive_binary_search(int arr[], int left, int right, int key) {
    if (right >= left) {
        int mid = left + (right-left)/2;
        if (arr[mid] == key) return mid;
        if (arr[mid] > key) return recursive_binary_search(arr, left, mid-1, key);
        return recursive_binary_search(arr, mid+1, right, key);
    }
    return -1;
}

/* -------------------- Number system functions -------------------- */
bool is_prime(int num) {
    if (num <= 1) return false;
    if (num <= 3) return true;
    if (num % 2 == 0) return false;
    for (int i = 3; i * i <= num; i += 2) {
        if (num % i == 0) return false;
    }
    return true;
}

long long factorial(int n) {
    if (n < 0) return -1; // indicate error
    long long res = 1;
    for (int i = 2; i <= n; ++i) res *= i;
    return res;
}

void generate_fibonacci(int n, int *fib) {
    if (n <= 0) return;
    fib[0] = 0;
    if (n == 1) return;
    fib[1] = 1;
    for (int i = 2; i < n; i++) fib[i] = fib[i-1] + fib[i-2];
}

int gcd(int a, int b) {
    if (a < 0) a = -a;
    if (b < 0) b = -b;
    while (b != 0) {
        int t = b;
        b = a % b;
        a = t;
    }
    return a;
}
int lcm(int a, int b) {
    if (a == 0 || b == 0) return 0;
    return (a / gcd(a,b)) * b;
}

/* -------------------- IPC demo (pipe, msg queue, shared memory + semaphore) -------------------- */
struct msg_buffer {
    long msg_type;
    char msg_text[100];
};

typedef struct {
    sem_t sem;        // semaphore for sync
    int fib[100];     // shared fibonacci
    int n;            // number of terms
} shm_segment_t;


/* -------------------- Device-driver simulation (updated) -------------------- */
#define RBUF_SIZE 8
typedef struct {
    char *buf[RBUF_SIZE];
    int head, tail;
    int count;
    pthread_mutex_t lock;
    pthread_cond_t not_empty;
    pthread_cond_t not_full;
} ring_buffer_t;

void ring_init(ring_buffer_t *r) {
    r->head = r->tail = r->count = 0;
    pthread_mutex_init(&r->lock, NULL);
    pthread_cond_init(&r->not_empty, NULL);
    pthread_cond_init(&r->not_full, NULL);
    for (int i=0;i<RBUF_SIZE;i++) r->buf[i] = NULL;
}

void ring_push(ring_buffer_t *r, char *s) {
    pthread_mutex_lock(&r->lock);
    while (r->count == RBUF_SIZE) pthread_cond_wait(&r->not_full, &r->lock);
    r->buf[r->tail] = strdup(s);
    r->tail = (r->tail + 1) % RBUF_SIZE;
    r->count++;
    pthread_cond_signal(&r->not_empty);
    pthread_mutex_unlock(&r->lock);
}

char* ring_pop(ring_buffer_t *r) {
    pthread_mutex_lock(&r->lock);
    while (r->count == 0) pthread_cond_wait(&r->not_empty, &r->lock);
    char *s = r->buf[r->head];
    r->buf[r->head] = NULL;
    r->head = (r->head + 1) % RBUF_SIZE;
    r->count--;
    pthread_cond_signal(&r->not_full);
    pthread_mutex_unlock(&r->lock);
    return s;
}

ring_buffer_t device_ring;
pid_t main_pid_for_signal = 0;

/* Optional SIGUSR1 handler (kept, but SIGUSR1 is optional in writer) */
void sigusr1_handler(int signo) {
    time_t now; time(&now);
    char *ts = ctime(&now);
    if (ts) ts[strcspn(ts, "\n")] = 0;
    printf("[DeviceDriver] Interrupt (SIGUSR1) at %s\n", ts ? ts : "unknown");
}

/* Writer thread */
void* device_writer_thread(void *arg) {
    int i = 1;
    pthread_setcanceltype(PTHREAD_CANCEL_DEFERRED, NULL);
    while (1) {
        char buff[128];
        snprintf(buff, sizeof(buff), "device-message-%d", i++);
        ring_push(&device_ring, buff);
       // printf("[DeviceDriver] Writer: Device producing data... ongoing simulation\n");
        /* If you want actual interrupt notifications, enable the kill below.
           By default left commented to avoid spamming SIGUSR1 during interactive input.
        */
        // if (main_pid_for_signal > 0) kill(main_pid_for_signal, SIGUSR1);
        sleep(1); /* cancellation point */
    }
    return NULL;
}

/* Reader thread */
void* device_reader_thread(void *arg) {
    pthread_setcanceltype(PTHREAD_CANCEL_DEFERRED, NULL);
    while (1) {
        char *msg = ring_pop(&device_ring); /* may block, cancellation point */
        if (msg) {
          //  printf("[DeviceDriver] Reader: Consumed %s\n", msg);
            free(msg);
        }
    }
    return NULL;
}

/* Friendly message for unimplemented ops */
void device_not_implemented(const char* feature) {
    printf("[DeviceDriver] Operation '%s' not implemented (requires real controller/processor)\n", feature);
}


/* -------------------- Bitwise operations -------------------- */
int set_bit(int num, int pos) {
    return num | (1 << pos);
}
int clear_bit(int num, int pos) {
    return num & ~(1 << pos);
}
int toggle_bit(int num, int pos) {
    return num ^ (1 << pos);
}
int check_bit(int num, int pos) {
    return (num >> pos) & 1;
}
int count_set_bits(int num) {
    int count = 0;
    while (num) {
        count += num & 1;
        num >>= 1;
    }
    return count;
}

/* -------------------- Data Structures Implementations -------------------- */
/* Queue (simple array-based) */
#define MAXQ 50
int qarr[MAXQ], qfront=-1, qrear=-1;
void enqueue(int val){
    if(qrear==MAXQ-1) { printf("Queue full!\n"); return; }
    if(qfront==-1) qfront=0;
    qarr[++qrear]=val;
}
void dequeue(){
    if(qfront==-1 || qfront>qrear){ printf("Queue empty!\n"); return; }
    printf("Dequeued: %d\n", qarr[qfront++]);
    if (qfront>qrear) { qfront = qrear = -1; }
}
void display_queue(){
    if(qfront==-1 || qfront>qrear){ printf("Empty queue\n"); return; }
    for(int i=qfront;i<=qrear;i++) printf("%d ",qarr[i]);
    printf("\n");
}
void queue_demo(){
    int ch,val;
    while(1){
        printf("\nQueue Menu: 1.Enqueue 2.Dequeue 3.Display 4.Back\n");
        scanf("%d",&ch);
        if(ch==1){ printf("Val: "); scanf("%d",&val); enqueue(val);}
        else if(ch==2) dequeue();
        else if(ch==3) display_queue();
        else break;
    }
}

/* Singly Linked List */
typedef struct Node{ int data; struct Node* next;}Node;
Node* head = NULL;
void insert_ll(int val){
    Node* n = malloc(sizeof(Node)); n->data = val; n->next = head; head = n;
}
void display_ll(){
    Node* t = head;
    while(t){ printf("%d->",t->data); t=t->next; }
    printf("NULL\n");
}
void delete_ll(int val){
    Node *t=head,*p=NULL;
    while(t && t->data!=val){ p=t; t=t->next;}
    if(!t){printf("Not found\n");return;}
    if(!p) head=head->next; else p->next=t->next;
    free(t);
}
void linkedlist_demo(){
    int ch,val;
    while(1){
        printf("\nLL Menu: 1.Insert 2.Delete 3.Display 4.Back\n");
        scanf("%d",&ch);
        if(ch==1){printf("Val: "); scanf("%d",&val); insert_ll(val);}
        else if(ch==2){printf("Val: "); scanf("%d",&val); delete_ll(val);}
        else if(ch==3) display_ll();
        else break;
    }
}

/* Doubly Linked List */
typedef struct DNode{ int data; struct DNode*prev,*next;}DNode;
DNode* dhead = NULL;
void insert_dl(int val){
    DNode* n = malloc(sizeof(DNode)); n->data = val; n->prev = NULL; n->next = dhead;
    if(dhead) dhead->prev = n; dhead = n;
}
void display_dl(){
    DNode* t = dhead; while(t){ printf("%d<->",t->data); t=t->next;} printf("NULL\n");
}
void delete_dl(int val){
    DNode*t=dhead; while(t && t->data!=val) t=t->next;
    if(!t){printf("Not found\n");return;}
    if(t->prev) t->prev->next = t->next; else dhead = t->next;
    if(t->next) t->next->prev = t->prev;
    free(t);
}
void doublylist_demo(){
    int ch,val;
    while(1){
        printf("\nDLL Menu: 1.Insert 2.Delete 3.Display 4.Back\n");
        scanf("%d",&ch);
        if(ch==1){printf("Val: "); scanf("%d",&val); insert_dl(val);}
        else if(ch==2){printf("Val: "); scanf("%d",&val); delete_dl(val);}
        else if(ch==3) display_dl();
        else break;
    }
}

/* Circular Linked List */
typedef struct CNode{int data; struct CNode* next;}CNode;
CNode* ctail = NULL;
void insert_cl(int val){
    CNode* n = malloc(sizeof(CNode)); n->data = val;
    if(!ctail){ ctail = n; n->next = n; }
    else { n->next = ctail->next; ctail->next = n; ctail = n; }
}
void display_cl(){
    if(!ctail){printf("Empty\n");return;}
    CNode* t = ctail->next;
    do{ printf("%d->",t->data); t = t->next; } while(t != ctail->next);
    printf("(back to head)\n");
}
void circularlist_demo(){
    int ch,val;
    while(1){
        printf("\nCLL Menu: 1.Insert 2.Display 3.Back\n");
        scanf("%d",&ch);
        if(ch==1){printf("Val: "); scanf("%d",&val); insert_cl(val);}
        else if(ch==2) display_cl();
        else break;
    }
}

/* Heap (max-heap) */
#define HMAX 100
int heap_arr[HMAX];
int hsize = 0;
void heap_insert(int val){
    if(hsize >= HMAX-1){ printf("Heap full\n"); return; }
    int i = ++hsize;
    heap_arr[i] = val;
    while (i > 1 && heap_arr[i/2] < heap_arr[i]) {
        int tmp = heap_arr[i]; heap_arr[i] = heap_arr[i/2]; heap_arr[i/2] = tmp;
        i /= 2;
    }
}
void display_heap(){
    if (hsize == 0) { printf("Heap empty\n"); return; }
    for(int i=1;i<=hsize;i++) printf("%d ", heap_arr[i]);
    printf("\n");
}
void heap_demo(){
    int ch,val;
    while(1){
        printf("\nHeap Menu: 1.Insert 2.Display 3.Back\n");
        scanf("%d",&ch);
        if(ch==1){printf("Val: "); scanf("%d",&val); heap_insert(val);}
        else if(ch==2) display_heap();
        else break;
    }
}

/* Binary Search Tree (simple) */
typedef struct TNode{int data; struct TNode*l,*r;}TNode;
TNode* root = NULL;
TNode* newNode(int v){ TNode* n = malloc(sizeof(TNode)); n->data = v; n->l = n->r = NULL; return n; }
TNode* insert_bst(TNode* r, int v) {
    if(!r) return newNode(v);
    if(v < r->data) r->l = insert_bst(r->l, v);
    else r->r = insert_bst(r->r, v);
    return r;
}
void inorder(TNode* r){ if(r){ inorder(r->l); printf("%d ", r->data); inorder(r->r);} }
void preorder(TNode* r){ if(r){ printf("%d ", r->data); preorder(r->l); preorder(r->r);} }
void postorder(TNode* r){ if(r){ postorder(r->l); postorder(r->r); printf("%d ", r->data);} }
void tree_demo(){
    int ch,val;
    while(1){
        printf("\nTree Menu: 1.Insert 2.Inorder 3.Preorder 4.Postorder 5.Back\n");
        scanf("%d",&ch);
        if(ch==1){printf("Val: "); scanf("%d",&val); root = insert_bst(root, val);}
        else if(ch==2){ inorder(root); printf("\n"); }
        else if(ch==3){ preorder(root); printf("\n"); }
        else if(ch==4){ postorder(root); printf("\n"); }
        else break;
    }
}

/* Graph (adjacency matrix, BFS/DFS) */
#define MAXV 20
int graph[MAXV][MAXV];
int V = 0;
void add_edge(int u, int v){ if(u<0||u>=V||v<0||v>=V) { printf("Bad vertices\n"); return; } graph[u][v] = 1; /* for undirected: */ graph[v][u] = 1; }
void bfs(int s){
    if(s<0||s>=V){ printf("Bad start\n"); return; }
    int visited[MAXV] = {0}, q[100], f=0, r=0;
    visited[s] = 1; q[r++] = s;
    while(f<r){
        int u = q[f++]; printf("%d ", u);
        for(int v=0; v<V; v++) if(graph[u][v] && !visited[v]) { visited[v]=1; q[r++]=v; }
    }
    printf("\n");
}
void dfs_util(int u, int visited[]){
    visited[u] = 1; printf("%d ", u);
    for(int v=0; v<V; v++) if(graph[u][v] && !visited[v]) dfs_util(v, visited);
}
void graph_demo(){
    int ch,u,v;
    printf("Enter number of vertices (<=%d): ", MAXV); scanf("%d",&V);
    if(V <= 0 || V > MAXV){ printf(BLINK BOLD_RED"Invalid V\n"RESET); return; }
    memset(graph, 0, sizeof(graph));
    while(1){
        printf("\nGraph Menu: 1.Add Edge 2.BFS 3.DFS 4.Back\n");
        scanf("%d",&ch);
        if(ch==1){ printf("u v: "); scanf("%d %d",&u,&v); add_edge(u,v); }
        else if(ch==2){ printf("Start: "); scanf("%d",&u); bfs(u); }
        else if(ch==3){ int visited[MAXV] = {0}; printf("Start: "); scanf("%d",&u); dfs_util(u, visited); printf("\n"); }
        else break;
    }
}

/* -------------------- Main Menu -------------------- */
int main()
{
    signal(SIGINT, interrupt_handler);  // Ctrl+C handler
    signal(SIGUSR1, sigusr1_handler);   // Device interrupt simulation

    main_pid_for_signal = getpid();
    ring_init(&device_ring);

    // start device simulation threads
    pthread_t writer_tid, reader_tid;
    pthread_create(&writer_tid, NULL, device_writer_thread, NULL);
    pthread_create(&reader_tid, NULL, device_reader_thread, NULL);

    int Number , option, n = 5 ;
    char choice;
    int arr[200];

    while(1){
		Device:
        printf("\n\n------------------ All In One C Demo ------------------\n");
        printf("Enter the Option for Perform the Operations\n");
        printf("1. Number_System\n2. Algorithm\n3. IPC\n4. DEVICE_DRIVER\n5. Polling Demo\n6. Bitwise Operations\n7. Data Structures\n8. EXIT\n");
        printf("Select: ");
        if (scanf("%d", &option) != 1) { // handle non-number input
            while (getchar() != '\n');
            printf(BLINK BOLD_RED"Invalid input\n"RESET);
            continue;
        }
        getchar(); // consume newline

        switch(option){

            case 1 : { // Number System
                while(1){
                    printf("\nNumber System Menu:\n");
                    printf("a. EVEN_ODD\nb. Factorial\nc. Prime Check\nd. Fibonacci\ne. GCD/LCM\nf. Back\n");
                    printf("choice: ");
                    if (scanf(" %c", &choice) != 1) { while(getchar()!='\n'); continue; }
                    getchar();

                    if (choice == 'a') {
                        printf("Enter the Number for checking EVEN|ODD: ");
                        scanf("%d", &Number); getchar();
                        if (Number == 0) printf("Number is Even: %d\n", Number);
                        else if (Number > 0) printf("Given Number is %s : %d\n", (Number & 1) ? "Odd" : "Even", Number);
                        else printf("Negative number provided: %d\n", Number);
                    } else if (choice == 'b') {
                        printf("Enter n for factorial (<=20 recommended): ");
                        scanf("%d", &Number); getchar();
                        long long f = factorial(Number);
                        if (f < 0) printf(BLINK BOLD_RED"Invalid (negative) input\n"RESET);
                        else printf("Factorial(%d) = %lld\n", Number, f);
                    } else if (choice == 'c') {
                        printf("Enter number to test prime: ");
                        scanf("%d", &Number); getchar();
                        printf("%d is %s\n", Number, is_prime(Number) ? "Prime" : "Not Prime");
                    } else if (choice == 'd') {
                        printf("Enter n (count) for Fibonacci (<=100): ");
                        scanf("%d", &Number); getchar();
                        if (Number > 100) Number = 100;
                        int fib[100] = {0};
                        generate_fibonacci(Number, fib);
                        printf("Fibonacci series: ");
                        for (int i=0;i<Number;i++) printf("%d ", fib[i]);
                        printf("\n");
                    } else if (choice == 'e') {
                        int a,b;
                        printf("Enter two numbers (a b): ");
                        scanf("%d %d", &a, &b); getchar();
                        printf("GCD(%d,%d)=%d, LCM=%d\n", a, b, gcd(a,b), lcm(a,b));
                    } else if (choice == 'f') {
						printf(BLINK BOLD_YELLOW"Back from the Number system ...\n"RESET);
                        break;
                    } else {
                        printf(BLINK BOLD_RED"Invalid option\n"RESET);
                    }
                }
                break;
            }

            case 2 : { // Algorithm
                while(1){
                    printf("\nAlgorithm Menu:\n");
                    printf("a. Searching_Algorithm\nb. Sorting_Algorithm\nc. Back\n");
                    printf("choice: ");
                    if (scanf(" %c", &choice) != 1) { while(getchar()!='\n'); continue; }
                    getchar();
                    if (choice == 'a') {
                        int select;
                        while(1){
                            printf("\nSearch Menu:\n1. Binary_Search\n2. Linear_search\n3. Recursive Binary Search\n4. Back\nSelect: ");
                            if (scanf("%d", &select) != 1) { while(getchar()!='\n'); continue; }
                            getchar();
                            if (select == 1) {
                                printf("Enter number of elements: ");
                                scanf("%d", &n); getchar();
                                inputArray(arr, n);
                                // binary search needs sorted array: ensure sorted
                                insertionSort(arr, n);
                                printf("Sorted array: "); printArray(arr,n);
                                int x;
                                printf("Enter element to search: ");
                                scanf("%d", &x); getchar();
                                int result = binary_search(arr, 0, n - 1, x);
                                if (result == -1) printf("Element not present\n");
                                else printf("Element present at index: %d\n", result);
                            } else if (select == 2) {
                                printf("Enter number of elements: ");
                                scanf("%d", &n); getchar();
                                inputArray(arr, n);
                                int X;
                                printf("Enter element to search: ");
                                scanf("%d", &X); getchar();
                                int res = linear_search(arr, n, X);
                                if (res == -1) printf("Element not present\n");
                                else printf("Element present at index: %d\n", res);
                            } else if (select == 3) {
                                printf("Enter number of elements: ");
                                scanf("%d", &n); getchar();
                                inputArray(arr, n);
                                insertionSort(arr, n);
                                printf("Sorted array: "); printArray(arr,n);
                                int key;
                                printf("Enter element to search: ");
                                scanf("%d",&key); getchar();
                                int idx = recursive_binary_search(arr, 0, n-1, key);
                                if (idx==-1) printf("Not found\n"); else printf("Found at %d\n", idx);
                            } else if (select == 4){
								printf(BLINK BOLD_YELLOW"Back from Seraching Algorithm Menu...\n "RESET);
								break;
                            } else printf(BLINK BOLD_RED"Invalid\n"RESET);
                        }
                    } else if (choice == 'b') {
                        while(1){
                            int sort_select;
                            printf("\nSort Menu:\n1.Selection 2.Bubble 3.Insertion 4.Merge 5.Quick 6.Heap 7.Counting 8.Shell 9.Radix 10.Cocktail 11.Bucket 12.Back\nSelect: ");
                            if (scanf("%d", &sort_select) != 1) { while(getchar()!='\n'); continue; }
                            getchar();
						   if (sort_select == 12) {
							   break;
						   } else if (sort_select >= 1 && sort_select <= 12) {
                                if (sort_select == 11) {
                                    int fn;
                                    printf("Enter number of float elements (<=50): ");
                                    scanf("%d",&fn); getchar();
                                    if (fn > 50) fn = 50;
                                    float farr[50];
                                    printf("Enter floats (0<=x<1 recommended):\n");
                                    for (int i=0;i<fn;i++) scanf("%f", &farr[i]);
                                    bucketSort(farr, fn);
                                    printf("Sorted floats: ");
                                    for (int i=0;i<fn;i++) printf("%f ", farr[i]);
                                    printf("\n");
                                } else {
                                    printf("Enter number of elements: ");
                                    scanf("%d", &n); getchar();
                                    inputArray(arr, n);
                                    printf("Original: "); printArray(arr,n);
                                    switch(sort_select) {
                                        case 1: selectionSort(arr,n); break;
                                        case 2: bubbleSort(arr,n); break;
                                        case 3: insertionSort(arr,n); break;
                                        case 4: mergeSort(arr,0,n-1); break;
                                        case 5: quickSort(arr,0,n-1); break;
                                        case 6: heapSort(arr,n); break;
                                        case 7: countSort(arr,n); break;
                                        case 8: shellSort(arr,n); break;
                                        case 9: radixSort(arr,n); break;
                                        case 10: cocktailSort(arr,n); break;
                                    }
                                    printf("Sorted: "); printArray(arr,n);
                                }
                            } 
                            else printf(BLINK BOLD_RED"Invalid\n"RESET);
                        }
                    } else if (choice == 'c'){
						printf(BLINK BOLD_YELLOW"Back from Sorting_Algorithm Menu ....\n"RESET);

						break;
                    } else printf(BLINK BOLD_RED"Invalid\n"RESET);
                }
                break;
            }

            case 3 : { // IPC
                printf("\n--- IPC Demo: Pipe + Message Queue + Shared Memory (with semaphore) ---\n");
                pid_t pid;
                int pipefd[2];
                char readbuffer[80];
                int n_terms, range;

                printf("Enter number of Fibonacci terms (<=100): ");
                scanf("%d", &n_terms);
                if (n_terms > 100) n_terms = 100;
                printf("Enter range for even/odd and prime check: ");
                scanf("%d", &range);

                if (pipe(pipefd) == -1) { perror("pipe"); break; }

                // SysV Message queue
                key_t key = ftok(".", 'A');
                int msgid = msgget(key, 0666 | IPC_CREAT);
                struct msg_buffer message;

                // Shared memory for struct (contains semaphore + fib[])
                int shmid = shmget(key, sizeof(shm_segment_t), 0666 | IPC_CREAT);
                if (shmid < 0) { perror("shmget"); break; }
                shm_segment_t *shm_ptr = (shm_segment_t*) shmat(shmid, NULL, 0);
                if (shm_ptr == (void*)-1) { perror("shmat"); break; }

                // Initialize semaphore in shared memory (only once from parent)
                if (sem_init(&shm_ptr->sem, 1, 0) != 0) {
                    perror("sem_init");
                }
                shm_ptr->n = n_terms;

                pid = fork();
                if (pid < 0) { perror("fork"); break; }

                if (pid > 0) { // Parent
                    close(pipefd[0]); // close read

                    // Write integers to pipe (1..range)
                    for (int i = 1; i <= range; i++) {
                        write(pipefd[1], &i, sizeof(int));
                    }
                    close(pipefd[1]);

                    // Send prime/non-prime via message queue
                    for (int i = 1; i <= range; i++) {
                        message.msg_type = 1;
                        if (is_prime(i)) sprintf(message.msg_text, "%d is Prime", i);
                        else sprintf(message.msg_text, "%d is Non-Prime", i);
                        msgsnd(msgid, &message, sizeof(message.msg_text), 0);
                    }

                    // Generate fibonacci into shared memory
                    generate_fibonacci(n_terms, shm_ptr->fib);

                    // Post semaphore to allow child to read Fibonacci
                    sem_post(&shm_ptr->sem);

                    wait(NULL);

                    // Cleanup
                    msgctl(msgid, IPC_RMID, NULL);
                    shmdt(shm_ptr);
                    shmctl(shmid, IPC_RMID, NULL);
                } else { // Child
                    close(pipefd[1]); // close write
                    // Read from pipe
                    printf("\n[Child] Even/Odd from pipe:\n");
                    for (int i = 1; i <= range; i++) {
                        int val;
                        if (read(pipefd[0], &val, sizeof(int)) > 0) {
                            if (val % 2 == 0) printf("%d is Even\n", val);
                            else printf("%d is Odd\n", val);
                        }
                    }
                    close(pipefd[0]);

                    // Read from message queue
                    printf("\n[Child] Prime checks from message queue:\n");
                    while (msgrcv(msgid, &message, sizeof(message.msg_text), 1, IPC_NOWAIT) != -1) {
                        printf("%s\n", message.msg_text);
                    }
                    // Wait for semaphore then read fibonacci from shared memory
                    printf("\n[Child] Waiting for shared memory (fibonacci)...\n");
                    sem_wait(&shm_ptr->sem);
                    printf("[Child] Fibonacci series received from shared memory: ");
                    for (int i = 0; i < shm_ptr->n; i++) printf("%d ", shm_ptr->fib[i]);
                    printf("\n");

                    shmdt(shm_ptr);
                    exit(0);
                }
                break;
            }

            case 4 : { /* DEVICE_DRIVER - start/stop inside menu */
                printf("\n[DeviceDriver] Simulation starting on request...\n");
                ring_init(&device_ring);
                main_pid_for_signal = getpid();
                signal(SIGUSR1, sigusr1_handler);

                pthread_t writer, reader;
                if (pthread_create(&writer, NULL, device_writer_thread, NULL) != 0) {
                    perror("pthread_create writer");
                    break;
                }
                if (pthread_create(&reader, NULL, device_reader_thread, NULL) != 0) {
                    perror("pthread_create reader");
                    pthread_cancel(writer);
                    pthread_join(writer, NULL);
                    break;
                }

                printf("[DeviceDriver] Threads started (Writer + Reader).\n");
                printf("[DeviceDriver] Press 'b' then Enter to stop simulation and return to main menu.\n");

                char ch;
                while (1) {
                    if (scanf(" %c", &ch) != 1) { while(getchar()!='\n'); continue; }
                    if (ch == 'b') {
                        printf("[DeviceDriver] Stopping simulation and returning to main menu...\n");
						printf(BLINK
						BOLD_PURPLE"Return to Main Menu from Device Driver...\n"RESET);
						goto Device;
                        pthread_cancel(writer);
                        pthread_cancel(reader);
                        pthread_join(writer, NULL);
                        pthread_join(reader, NULL);
                        /* drain ring buffer to free memory if desired */
                        while (device_ring.count > 0) {
                            char *m = ring_pop(&device_ring);
                            if (m) free(m);
                        }
                        break;
                    } else {
                        printf(BLINK BOLD_RED"[DeviceDriver] Unknown/Invalid input '%c'. Return to Main Menu...\n"RESET, ch);
						goto Device;
                    }
                }
                break;
            }

            case 6: { // Bitwise Operations
                while (1) {
                    printf("\nBitwise Operations Menu:\n");
                    printf("a. Set a bit\nb. Clear a bit\nc. Toggle a bit\nd. Check a bit\ne. Count set bits\nf. Back\n");
                    printf("choice: ");
                    if (scanf(" %c", &choice) != 1) { while(getchar()!='\n'); continue; }
                    getchar();

                    int num, pos;
                    switch (choice) {
                        case 'a':
                            printf("Enter number and bit position: ");
                            scanf("%d %d", &num, &pos); getchar();
                            printf("Result = %d\n", set_bit(num,pos));
                            break;
                        case 'b':
                            printf("Enter number and bit position: ");
                            scanf("%d %d", &num, &pos); getchar();
                            printf("Result = %d\n", clear_bit(num,pos));
                            break;
                        case 'c':
                            printf("Enter number and bit position: ");
                            scanf("%d %d", &num, &pos); getchar();
                            printf("Result = %d\n", toggle_bit(num,pos));
                            break;
                        case 'd':
                            printf("Enter number and bit position: ");
                            scanf("%d %d", &num, &pos); getchar();
                            printf("Bit %d is %s\n", pos, check_bit(num,pos) ? "SET" : "CLEAR");
                            break;
                        case 'e':
                            printf("Enter number: ");
                            scanf("%d", &num); getchar();
                            printf("Set bits count = %d\n", count_set_bits(num));
                            break;
                        case 'f':
                            goto back_bitwise;
                        default:
                            printf(BLINK BOLD_RED"Invalid option\n"RESET);
                    }
                }
                back_bitwise:
                break;
            }

            case 7: { // Data Structures
                while (1) {
                    printf("\nData Structures Menu:\n");
                    printf("a. Queue\nb. Singly Linked List\nc. Doubly Linked List\nd. Circular Linked List\ne. Heap\nf. Binary Tree\ng. Graph\nh. Back\n");
                    printf("choice: ");
                    if (scanf(" %c", &choice) != 1) { while(getchar()!='\n'); continue; }
                    getchar();

                    switch(choice) {
                        case 'a': queue_demo(); break;
                        case 'b': linkedlist_demo(); break;
                        case 'c': doublylist_demo(); break;
                        case 'd': circularlist_demo(); break;
                        case 'e': heap_demo(); break;
                        case 'f': tree_demo(); break;
                        case 'g': graph_demo(); break;
                        case 'h': goto back_ds;
                        default: printf(BLINK BOLD_RED"Invalid\n"RESET);
                    }
                }
                back_ds:
                break;
            }

            case 8:
                printf(BLINK BOLD_GREEN"Exiting...\n"RESET);
                exit(0);

            default:
                printf(BLINK BOLD_RED"Invalid Option\n"RESET);
                break;
        }
    }

    return 0;
}

