#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <cmath>
#include <mpi.h>
#include <stddef.h>
#include <string.h>
#include <float.h>
#include <time.h>

struct Point{
    float coord[2];
    int index;
};

struct PointDomain{
    int i;
    int j;
    float coord[2];
    int domain;
};

struct Comparator{
    int first;
    int second;
    Comparator(int f, int s): first(f), second(s) {}
};

int axis = 1;
std::vector<Comparator> comparators;
MPI_Datatype MPI_POINT_TYPE;
MPI_Datatype MPI_POINT_DOMAIN_TYPE;

int compare_points(const void *a, const void *b)
{
    if(((Point*)a)->coord[axis] < ((Point*)b)->coord[axis]){
        return -1;
    }else if(((Point*)a)->coord[axis] > ((Point*)b)->coord[axis]){
        return 1;
    }else{
        return 0;
    }
}

void join(int *array1, int count1, int *array2, int count2)
{
    if(count1 + count2 == 1){
        return;
    }else if(count1 + count2 == 2){
        if(count1 == 2)
            comparators.push_back(Comparator(array1[0], array1[1]));
        else if(count2 == 2)
            comparators.push_back(Comparator(array2[0], array2[1]));
        else
            comparators.push_back(Comparator(array1[0], array2[0]));
    }else{
        int oddcount1 = count1 / 2;
        int evencount1 = count1 - oddcount1;
        int *even1 = new int[evencount1];
        int *odd1 = new int[oddcount1];
        for(int i = 0, j = 0, k = 0; i < count1; i++){
            if(i % 2 == 0){
                even1[j] = array1[i];
                j++;
            }else{
                odd1[k] = array1[i];
                k++;
            }
        }

        int oddcount2 = count2 / 2;
        int evencount2 = count2 - oddcount2;
        int *even2 = new int[evencount2];
        int *odd2 = new int[oddcount2];
        for(int i = 0, j = 0, k = 0; i < count2; i++){
            if(i % 2 == 0){
                even2[j] = array2[i];
                j++;
            }else{
                odd2[k] = array2[i];
                k++;
            }
        }

        join(even1, evencount1, even2, evencount2);
        join(odd1, oddcount1, odd2, oddcount2);
        int *res = new int[count1 + count2];
        for(int i = 0; i < count1; i++)
            res[i] = array1[i];
        for(int i = 0; i < count2; i++)
            res[count1 + i] = array2[i];
        for(int i = 1; i < count1 + count2 - 1; i += 2)
            comparators.push_back(Comparator(res[i], res[i+1]));
        delete [] res;
        delete [] odd1;
        delete [] even1;
        delete [] odd2;
        delete [] even2;
    }
}

void sort(int *array, int count)
{
    int size1 = count / 2;
    int size2 = count - size1;
    if(count > 1){
        int *array1 = new int[size1];
        int *array2 = new int[size2];
        for(int i = 0; i < size1; i++)
            array1[i] = array[i];
        for(int i = 0; i < size2; i++)
            array2[i] = array[size1 + i];
        sort(array1, size1);
        sort(array2, size2);
        join(array1, size1, array2, size2);
        delete [] array1;
        delete [] array2;
    }
}

void create_shedule(int start, int count)
{
    int *array = new int[count];
    for(int i = 0; i < count; i++)
        array[i] = start + i;
    sort(array, count);
    delete [] array;
}

float x(int i, int j)
{
    return ((float)rand()/(float)RAND_MAX * 100000.0 - 50000.0);
}

float y(int i, int j)
{
    return ((float)rand()/(float)RAND_MAX * 100000.0 - 50000.0);
}

void create_point_type()
{
    const int nitems = 2;
    int blocklengths[2] = {2, 1};
    MPI_Datatype types[2] = {MPI_FLOAT, MPI_INT};
    MPI_Aint offsets[2];

    offsets[0] = offsetof(Point, coord);
    offsets[1] = offsetof(Point, index);

    MPI_Type_create_struct(nitems, blocklengths, offsets, types,
            &MPI_POINT_TYPE);
    MPI_Type_commit(&MPI_POINT_TYPE);
}

void create_point_domain_type()
{
    const int nitems = 4;
    int blocklengths[4] = {1, 1, 2, 1};
    MPI_Datatype types[4] = {MPI_INT, MPI_INT, MPI_FLOAT, MPI_INT};
    MPI_Aint offsets[4];

    offsets[0] = offsetof(PointDomain, i);
    offsets[1] = offsetof(PointDomain, j);
    offsets[2] = offsetof(PointDomain, coord);
    offsets[3] = offsetof(PointDomain, domain);

    MPI_Type_create_struct(nitems, blocklengths, offsets, types,
            &MPI_POINT_DOMAIN_TYPE);
    MPI_Type_commit(&MPI_POINT_DOMAIN_TYPE);
}

Point* init_array(int n1, int n2, int rank, int proc_count)
{
    int elem_total = n1 * n2;
    int elem_per_proc_count = ceil(elem_total / (double)proc_count);
    int proc_without_fictive_elements = elem_total % proc_count;
    int temp = elem_total / proc_count;
    int elem_non_fictive = rank < proc_without_fictive_elements ?
        temp + 1 : temp;
    int start = rank < proc_without_fictive_elements ?
        rank * elem_per_proc_count : proc_without_fictive_elements *
        elem_per_proc_count + (rank - proc_without_fictive_elements) * temp;

    Point *array = new Point[elem_per_proc_count];
    int k = 0;
    for(; k < elem_non_fictive; k++){
        int i = (k + start) / n2;
        int j = (k + start) % n2;
        array[k].index = k + start;
        array[k].coord[0] = x(i, j);
        array[k].coord[1] = y(i, j);
    }
    for(; k < elem_per_proc_count; k++){
        array[k].index = -1;
        array[k].coord[0] = FLT_MAX;
        array[k].coord[1] = FLT_MAX;
    }
    return array;
}

void write_array(PointDomain *array, int n, char *filename, int n1, int n2,
        int rank)
{
    MPI_Status s;

    MPI_File output;
    if(MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_CREATE |
                MPI_MODE_WRONLY, MPI_INFO_NULL, &output) != MPI_SUCCESS){
        if(rank == 0)
            printf("File '%s' can't be opened.\n", filename);
        MPI_Finalize();
        exit(1);
    }
    MPI_File_set_size(output, 0);

    MPI_File_write_ordered(output, &n1, rank == 0, MPI_INT, &s);
    MPI_File_write_ordered(output, &n2, rank == 0, MPI_INT, &s);
    MPI_File_write_ordered(output, array, n, MPI_POINT_DOMAIN_TYPE, &s);
    
    MPI_File_close(&output);
}

void sort_array(Point*& array, int elem_per_proc_count, MPI_Comm communicator)
{
    MPI_Status s;
    int rank;
    MPI_Comm_rank(communicator, &rank);
    Point *from_another_proc = new Point[elem_per_proc_count];
    Point *tmp = new Point[elem_per_proc_count];
    qsort(array, elem_per_proc_count, sizeof(Point), compare_points);
    for(unsigned int ui = 0; ui < comparators.size(); ui++){
        Comparator comparator = comparators[ui];
        if(rank == comparator.first){
            MPI_Send(array, elem_per_proc_count, MPI_POINT_TYPE,
                    comparator.second, 0, communicator);
            MPI_Recv(from_another_proc, elem_per_proc_count, MPI_POINT_TYPE,
                    comparator.second, 0, communicator, &s);
            
            for(int i = 0, j = 0, k = 0; i < elem_per_proc_count; i++){
                Point a = array[j];
                Point b = from_another_proc[k];
                // a < b
                if(compare_points(&a, &b) == -1){
                    tmp[i] = a;
                    j++;
                }else{
                    tmp[i] = b;
                    k++;
                }
            }
            
            Point *temp = array;
            array = tmp;
            tmp = temp;
        }else if(rank == comparator.second){
            MPI_Recv(from_another_proc, elem_per_proc_count, MPI_POINT_TYPE,
                    comparator.first, 0, communicator, &s);
            MPI_Send(array, elem_per_proc_count, MPI_POINT_TYPE,
                    comparator.first, 0, communicator);
            
            int start = elem_per_proc_count - 1;
            for(int i = start, j = start, k = start; i >= 0; i--){
                Point a = array[j];
                Point b = from_another_proc[k];
                // a > b
                if(compare_points(&a, &b) == 1){
                    tmp[i] = a;
                    j--;
                }else{
                    tmp[i] = b;
                    k--;
                }
            }
            
            Point *temp = array;
            array = tmp;
            tmp = temp;
        }
    }
    delete [] from_another_proc;
    delete [] tmp;
}

void local_decompose(Point *a, int *domains, int domain_start, int k,
        int array_start, int n)
{
    if(k == 1){
        for(int i = 0; i < n; i++)
            domains[array_start + i] = domain_start;
        return;
    }

    axis = !axis;
    qsort(a + array_start, n, sizeof(Point), compare_points);

    int k1 = (k + 1) / 2;
    int k2 = k - k1;
    int n1 = n * (k1 / (double)k);
    int n2 = n - n1;
    local_decompose(a, domains, domain_start, k1, array_start, n1);
    local_decompose(a, domains, domain_start + k1, k2, array_start + n1, n2);
}

int remove_fictive(Point **array, int size)
{
    Point *res = new Point[size];
    int j = 0;
    for(int i = 0; i < size; i++){
        if((*array)[i].index == -1)
            continue;
        res[j++] = (*array)[i];
    }
    delete [] (*array);
    *array = res;
    return j;
}

Point* align(Point *array, int size, int proc_count)
{
    int epp = ceil(size / (double)proc_count);
    Point *res = new Point[epp * proc_count];
    int k = 0;
    for(int i = 0; i < size; i++, k++){
        res[k] = array[i];
        array[i].index = -1;
        array[i].coord[0] = FLT_MAX;
        array[i].coord[1] = FLT_MAX;
    }
    for(; k < epp * proc_count; k++){
        res[k].index = -1;
        res[k].coord[0] = FLT_MAX;
        res[k].coord[1] = FLT_MAX;
    }
    return res;
}

void decompose(Point **array, int **domains, int domain_start, int k, int n,
        int *elem_per_proc, MPI_Comm communicator)
{
    int rank, proc_count, actual_size;
    MPI_Comm_rank(communicator, &rank);
    MPI_Comm_size(communicator, &proc_count);
    if(proc_count == 1){
        actual_size = remove_fictive(array, *elem_per_proc);
        *domains = new int[actual_size];
        local_decompose(*array, *domains, domain_start, k, 0, actual_size);
        *elem_per_proc = actual_size;
        return;
    }
    if(k == 1){
        actual_size = remove_fictive(array, *elem_per_proc);
        *domains = new int[actual_size];
        for(int i = 0; i < actual_size; i++)
            (*domains)[i] = domain_start;
        *elem_per_proc = actual_size;
        return;
    }

    MPI_Status s;
    int k1 = (k + 1) / 2;
    int k2 = k - k1;
    int n1 = n * (k1 / (double)k);
    int n2 = n - n1;
    int pc = n1 / (*elem_per_proc);
    int middle = n1 % (*elem_per_proc);
    int color;
    if(pc == 0)
        color = rank > pc ? 0 : 1;
    else
        color = rank >= pc ? 0 : 1;

    axis = !axis;
    comparators.clear();
    create_shedule(0, proc_count);
    sort_array(*array, *elem_per_proc, communicator);
    MPI_Comm new_comm;
    MPI_Comm_split(communicator, color, rank, &new_comm);

    if(pc == 0){
        int epp = ceil(((*elem_per_proc) - middle) / (double)(proc_count - pc));
        if(rank == pc){
            Point *temp = align((*array) + middle, (*elem_per_proc) - middle,
                    proc_count - pc - 1);
            for(int i = pc + 1, j = 0; i < proc_count; i++, j++)
                MPI_Send(temp + j*epp, epp, MPI_POINT_TYPE, i, 0, communicator);
            delete [] temp;
            decompose(array, domains, domain_start, k1, n1, elem_per_proc,
                    new_comm);
        }else{
            Point *arr = new Point[*elem_per_proc + epp];
            MPI_Recv(arr + (*elem_per_proc), epp, MPI_POINT_TYPE, pc, 0,
                    communicator, &s);
            memcpy(arr, *array, (*elem_per_proc) * sizeof(Point));
            *elem_per_proc += epp;
            delete [] (*array);
            *array = arr;
            decompose(array, domains, domain_start + k1, k2, n2, elem_per_proc,
                    new_comm);
        }
        return;
    }
    if(rank <= pc){
        int epp = ceil(middle / (double)pc);
        if(rank == pc){
            Point *temp = align(*array, middle, pc);
            for(int i = 0; i < pc; i++)
                MPI_Send(temp + i*epp, epp, MPI_POINT_TYPE, i, 0, communicator);
            delete [] temp;
        }else{
            Point *arr = new Point[*elem_per_proc + epp];
            MPI_Recv(arr + (*elem_per_proc), epp, MPI_POINT_TYPE, pc, 0,
                    communicator, &s);
            memcpy(arr, *array, (*elem_per_proc) * sizeof(Point));
            *elem_per_proc += epp;
            delete [] (*array);
            *array = arr;
        }
    }
    if(rank < pc){
        decompose(array, domains, domain_start, k1, n1, elem_per_proc,
                new_comm);
    }else{
        decompose(array, domains, domain_start + k1, k2, n2, elem_per_proc,
                new_comm);
    }
}

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    int rank, proc_count;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &proc_count);
    if(argc != 5){
        if(rank == 0){
            printf("Usage: mpirun -np <proc_num> decompose <k> <n1> <n2> "
                    "<output file>\n");
        }
        MPI_Finalize();
        return 0;
    }
    create_point_type();
    create_point_domain_type();
    int k = strtol(argv[1], NULL, 10);
    int n1 = strtol(argv[2], NULL, 10);
    int n2 = strtol(argv[3], NULL, 10);
    int elem_per_proc = ceil(n1 * n2 / (double)proc_count);
    srand(time(NULL) + rank * proc_count);

    // initialization
    Point *array = init_array(n1, n2, rank, proc_count);
    int *domains;

    // decomposition
    decompose(&array, &domains, 0, k, n1 * n2, &elem_per_proc, MPI_COMM_WORLD);

    // writing
    int wcount = 0;
    PointDomain *res = new PointDomain[elem_per_proc];
    for(int i = 0; i < elem_per_proc; i++){
        if(array[i].index == -1)
            continue;
        res[wcount].coord[0] = array[i].coord[0];
        res[wcount].coord[1] = array[i].coord[1];
        res[wcount].i = array[i].index / n2;
        res[wcount].j = array[i].index % n2;
        res[wcount].domain = domains[i];
        wcount++;
    }
    delete [] domains;
    delete [] array;
    write_array(res, wcount, argv[4], n1, n2, rank);
    delete [] res;

    MPI_Finalize();
    return 0;
}
