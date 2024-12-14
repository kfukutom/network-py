#include <Python.h>
#include <structmember.h>
#include <math.h>
#include <float.h>

#define MAX_NODES 1000000
#define MAX_EDGES 5000000
#define EARTH_RADIUS 6371000.0 // in meters

typedef struct {
    double lat;
    double lon;
    int id;
} Node;

typedef struct {
    int source;
    int target;
    double length;
    int id;
} Edge;

typedef struct {
    PyObject_HEAD
    Node* nodes;
    int node_count;
    Edge* edges;
    int edge_count;
} UrbanNetwork;

static void
UrbanNetwork_dealloc(UrbanNetwork* self)
{
    free(self->nodes);
    free(self->edges);
    Py_TYPE(self)->tp_free((PyObject*)self);
}

static PyObject*
UrbanNetwork_new(PyTypeObject* type, PyObject* args, PyObject* kwds)
{
    UrbanNetwork* self;
    self = (UrbanNetwork*)type->tp_alloc(type, 0);
    if (self != NULL) {
        self->nodes = (Node*)malloc(MAX_NODES * sizeof(Node));
        self->edges = (Edge*)malloc(MAX_EDGES * sizeof(Edge));
        self->node_count = 0;
        self->edge_count = 0;
    }
    return (PyObject*)self;
}

static int
UrbanNetwork_init(UrbanNetwork* self, PyObject* args, PyObject* kwds)
{
    return 0;
}

static PyObject*
UrbanNetwork_add_node(UrbanNetwork* self, PyObject* args)
{
    double lat, lon;
    int id;
    if (!PyArg_ParseTuple(args, "ddi", &lat, &lon, &id)) {
        return NULL;
    }
    if (self->node_count >= MAX_NODES) {
        PyErr_SetString(PyExc_OverflowError, "Maximum number of nodes reached");
        return NULL;
    }
    self->nodes[self->node_count].lat = lat;
    self->nodes[self->node_count].lon = lon;
    self->nodes[self->node_count].id = id;
    self->node_count++;
    Py_RETURN_NONE;
}

static PyObject*
UrbanNetwork_add_edge(UrbanNetwork* self, PyObject* args)
{
    int source, target, id;
    double length;
    if (!PyArg_ParseTuple(args, "iidi", &source, &target, &length, &id)) {
        return NULL;
    }
    if (self->edge_count >= MAX_EDGES) {
        PyErr_SetString(PyExc_OverflowError, "Maximum number of edges reached");
        return NULL;
    }
    self->edges[self->edge_count].source = source;
    self->edges[self->edge_count].target = target;
    self->edges[self->edge_count].length = length;
    self->edges[self->edge_count].id = id;
    self->edge_count++;
    Py_RETURN_NONE;
}

static double
haversine_distance(double lat1, double lon1, double lat2, double lon2)
{
    double dlat = (lat2 - lat1) * M_PI / 180.0;
    double dlon = (lon2 - lon1) * M_PI / 180.0;
    double a = sin(dlat/2) * sin(dlat/2) + cos(lat1 * M_PI / 180.0) * cos(lat2 * M_PI / 180.0) * sin(dlon/2) * sin(dlon/2);
    double c = 2 * atan2(sqrt(a), sqrt(1-a));
    return EARTH_RADIUS * c;
}

static PyObject*
UrbanNetwork_calculate_distances(UrbanNetwork* self)
{
    for (int i = 0; i < self->edge_count; i++) {
        int source = self->edges[i].source;
        int target = self->edges[i].target;
        double distance = haversine_distance(
            self->nodes[source].lat, self->nodes[source].lon,
            self->nodes[target].lat, self->nodes[target].lon
        );
        self->edges[i].length = distance;
    }
    Py_RETURN_NONE;
}

static PyObject*
UrbanNetwork_get_node_degree(UrbanNetwork* self, PyObject* args)
{
    int node_id;
    if (!PyArg_ParseTuple(args, "i", &node_id)) {
        return NULL;
    }
    int degree = 0;
    for (int i = 0; i < self->edge_count; i++) {
        if (self->edges[i].source == node_id || self->edges[i].target == node_id) {
            degree++;
        }
    }
    return PyLong_FromLong(degree);
}

static PyObject*
UrbanNetwork_get_average_degree(UrbanNetwork* self)
{
    double total_degree = 0;
    for (int i = 0; i < self->node_count; i++) {
        int degree = 0;
        for (int j = 0; j < self->edge_count; j++) {
            if (self->edges[j].source == self->nodes[i].id || self->edges[j].target == self->nodes[i].id) {
                degree++;
            }
        }
        total_degree += degree;
    }
    return PyFloat_FromDouble(total_degree / self->node_count);
}

static PyObject*
UrbanNetwork_get_network_density(UrbanNetwork* self)
{
    double max_edges = self->node_count * (self->node_count - 1) / 2.0;
    return PyFloat_FromDouble(self->edge_count / max_edges);
}

static PyObject*
UrbanNetwork_get_average_path_length(UrbanNetwork* self)
{
    // Floyd-Warshall algorithm for average path length
    double** dist = (double**)malloc(self->node_count * sizeof(double*));
    for (int i = 0; i < self->node_count; i++) {
        dist[i] = (double*)malloc(self->node_count * sizeof(double));
        for (int j = 0; j < self->node_count; j++) {
            if (i == j) {
                dist[i][j] = 0;
            } else {
                dist[i][j] = DBL_MAX;
            }
        }
    }

    for (int i = 0; i < self->edge_count; i++) {
        int source = self->edges[i].source;
        int target = self->edges[i].target;
        dist[source][target] = self->edges[i].length;
        dist[target][source] = self->edges[i].length;
    }

    for (int k = 0; k < self->node_count; k++) {
        for (int i = 0; i < self->node_count; i++) {
            for (int j = 0; j < self->node_count; j++) {
                if (dist[i][k] + dist[k][j] < dist[i][j]) {
                    dist[i][j] = dist[i][k] + dist[k][j];
                }
            }
        }
    }

    double total_path_length = 0;
    int path_count = 0;
    for (int i = 0; i < self->node_count; i++) {
        for (int j = i + 1; j < self->node_count; j++) {
            if (dist[i][j] != DBL_MAX) {
                total_path_length += dist[i][j];
                path_count++;
            }
        }
    }

    for (int i = 0; i < self->node_count; i++) {
        free(dist[i]);
    }
    free(dist);

    return PyFloat_FromDouble(total_path_length / path_count);
}

static PyObject*
UrbanNetwork_get_betweenness_centrality(UrbanNetwork* self)
{
    // Brandes' algorithm for betweenness centrality
    double* centrality = (double*)calloc(self->node_count, sizeof(double));
    
    for (int s = 0; s < self->node_count; s++) {
        int* stack = (int*)malloc(self->node_count * sizeof(int));
        int stack_top = 0;
        
        double* dist = (double*)malloc(self->node_count * sizeof(double));
        int* pred = (int*)malloc(self->node_count * sizeof(int));
        double* sigma = (double*)malloc(self->node_count * sizeof(double));
        double* delta = (double*)calloc(self->node_count, sizeof(double));
        
        for (int t = 0; t < self->node_count; t++) {
            dist[t] = -1;
            pred[t] = -1;
            sigma[t] = 0;
        }
        
        dist[s] = 0;
        sigma[s] = 1;
        
        int* queue = (int*)malloc(self->node_count * sizeof(int));
        int queue_start = 0, queue_end = 0;
        queue[queue_end++] = s;
        
        while (queue_start != queue_end) {
            int v = queue[queue_start++];
            stack[stack_top++] = v;
            
            for (int i = 0; i < self->edge_count; i++) {
                if (self->edges[i].source == v || self->edges[i].target == v) {
                    int w = (self->edges[i].source == v) ? self->edges[i].target : self->edges[i].source;
                    if (dist[w] < 0) {
                        queue[queue_end++] = w;
                        dist[w] = dist[v] + 1;
                    }
                    if (dist[w] == dist[v] + 1) {
                        sigma[w] += sigma[v];
                        pred[w] = v;
                    }
                }
            }
        }
        
        while (--stack_top >= 0) {
            int w = stack[stack_top];
            if (pred[w] != -1) {
                delta[pred[w]] += (1 + delta[w]) * (sigma[pred[w]] / sigma[w]);
            }
            if (w != s) {
                centrality[w] += delta[w];
            }
        }
        
        free(stack);
        free(dist);
        free(pred);
        free(sigma);
        free(delta);
        free(queue);
    }
    
    PyObject* result = PyList_New(self->node_count);
    for (int i = 0; i < self->node_count; i++) {
        PyList_SetItem(result, i, PyFloat_FromDouble(centrality[i] / ((self->node_count - 1) * (self->node_count - 2))));
    }
    
    free(centrality);
    return result;
}

static PyObject*
UrbanNetwork_get_closeness_centrality(UrbanNetwork* self)
{
    PyObject* result = PyList_New(self->node_count);
    
    for (int i = 0; i < self->node_count; i++) {
        double* dist = (double*)malloc(self->node_count * sizeof(double));
        for (int j = 0; j < self->node_count; j++) {
            dist[j] = (i == j) ? 0 : DBL_MAX;
        }
        
        // Dijkstra's algorithm
        int* visited = (int*)calloc(self->node_count, sizeof(int));
        for (int count = 0; count < self->node_count - 1; count++) {
            int u = -1;
            double min_dist = DBL_MAX;
            for (int v = 0; v < self->node_count; v++) {
                if (!visited[v] && dist[v] < min_dist) {
                    min_dist = dist[v];
                    u = v;
                }
            }
            
            if (u == -1) break;
            
            visited[u] = 1;
            
            for (int j = 0; j < self->edge_count; j++) {
                if (self->edges[j].source == u || self->edges[j].target == u) {
                    int v = (self->edges[j].source == u) ? self->edges[j].target : self->edges[j].source;
                    if (!visited[v] && dist[u] + self->edges[j].length < dist[v]) {
                        dist[v] = dist[u] + self->edges[j].length;
                    }
                }
            }
        }
        
        double sum_dist = 0;
        int reachable_nodes = 0;
        for (int j = 0; j < self->node_count; j++) {
            if (dist[j] != DBL_MAX) {
                sum_dist += dist[j];
                reachable_nodes++;
            }
        }
        
        double closeness = (reachable_nodes > 1) ? (reachable_nodes - 1) / sum_dist : 0;
        PyList_SetItem(result, i, PyFloat_FromDouble(closeness));
        
        free(dist);
        free(visited);
    }
    
    return result;
}

static PyMethodDef UrbanNetwork_methods[] = {
    {"add_node", (PyCFunction)UrbanNetwork_add_node, METH_VARARGS, "Add a node to the network"},
    {"add_edge", (PyCFunction)UrbanNetwork_add_edge, METH_VARARGS, "Add an edge to the network"},
    {"calculate_distances", (PyCFunction)UrbanNetwork_calculate_distances, METH_NOARGS, "Calculate distances for all edges"},
    {"get_node_degree", (PyCFunction)UrbanNetwork_get_node_degree, METH_VARARGS, "Get the degree of a node"},
    {"get_average_degree", (PyCFunction)UrbanNetwork_get_average_degree, METH_NOARGS, "Get the average degree of the network"},
    {"get_network_density", (PyCFunction)UrbanNetwork_get_network_density, METH_NOARGS, "Get the density of the network"},
    {"get_average_path_length", (PyCFunction)UrbanNetwork_get_average_path_length, METH_NOARGS, "Get the average path length of the network"},
    {"get_betweenness_centrality", (PyCFunction)UrbanNetwork_get_betweenness_centrality, METH_NOARGS, "Get betweenness centrality for all nodes"},
    {"get_closeness_centrality", (PyCFunction)UrbanNetwork_get_closeness_centrality, METH_NOARGS, "Get closeness centrality for all nodes"},
    {NULL}  /* Sentinel */
};

static PyTypeObject UrbanNetworkType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    .tp_name = "urban_mapping.UrbanNetwork",
    .tp_doc = "Urban network object",
    .tp_basicsize = sizeof(UrbanNetwork),
    .tp_itemsize = 0,
    .tp_flags = Py_TPFLAGS_DEFAULT,
    .tp_new = UrbanNetwork_new,
    .tp_init = (initproc)UrbanNetwork_init,
    .tp_dealloc = (destructor)UrbanNetwork_dealloc,
    .tp_methods = UrbanNetwork_methods,
};