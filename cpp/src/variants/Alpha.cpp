#include "../../include/petls_headers/variants/Alpha.hpp"
#include <vector>
#include <set>
#include <iostream>
#ifdef PETLS_USE_ALPHA_COMPLEX



// Imports and typedefs for Gudhi
#include <gudhi/Simplex_tree.h>
#include <gudhi/Points_3D_off_io.h>
#include <boost/variant.hpp>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Alpha_shape_3.h>
#include <CGAL/Alpha_shape_vertex_base_3.h>
#include <CGAL/Alpha_shape_cell_base_3.h>
#include <CGAL/iterator.h>

// Alpha_shape_3 templates type definitions
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Alpha_shape_vertex_base_3<Kernel> Vb;
typedef CGAL::Alpha_shape_cell_base_3<Kernel> Fb;
typedef CGAL::Triangulation_data_structure_3<Vb, Fb> Tds;
typedef CGAL::Delaunay_triangulation_3<Kernel, Tds> Triangulation_3;
typedef CGAL::Alpha_shape_3<Triangulation_3> Alpha_shape_3;
// From file type definition
typedef Kernel::Point_3 Point;
// filtration with alpha values needed type definition
typedef Alpha_shape_3::FT Alpha_value_type;
typedef CGAL::Object Object;
typedef CGAL::Dispatch_output_iterator<
    CGAL::cpp11::tuple<Object, Alpha_value_type>,
    CGAL::cpp11::tuple<std::back_insert_iterator<std::vector<Object>>,
                       std::back_insert_iterator<std::vector<Alpha_value_type>>>>
    Dispatch;
typedef Alpha_shape_3::Cell_handle Cell_handle;
typedef Alpha_shape_3::Facet Facet;
typedef Alpha_shape_3::Edge Edge;
typedef std::list<Alpha_shape_3::Vertex_handle> Vertex_list;
// gudhi type definition
typedef Gudhi::Simplex_tree<> Simplex_tree;
typedef Simplex_tree::Vertex_handle Simplex_tree_vertex;
typedef std::map<Alpha_shape_3::Vertex_handle, Simplex_tree_vertex> Alpha_shape_simplex_tree_map;
typedef std::pair<Alpha_shape_3::Vertex_handle, Simplex_tree_vertex> Alpha_shape_simplex_tree_pair;
typedef std::vector<Simplex_tree_vertex> Simplex_tree_vector_vertex;

namespace petls{
    void alpha_points(std::vector<Point> &points, int dim_max, std::vector<std::vector<std::tuple<int,int,int>>>& boundaries_triples,std::vector<std::vector<filtration_type>>& filtrations);
    void alpha_OFF(const char* filename, int dim_max, std::vector<std::vector<std::tuple<int,int,int>>>& boundaries_triples,std::vector<std::vector<filtration_type>>& filtrations);
            

    Alpha::Alpha(const char* filename, int max_dim) : petls::Complex() { // call default constructor of parent class

        std::vector<std::vector<std::tuple<int,int,int>>> boundaries_triples(max_dim);
        std::vector<std::vector<filtration_type>> filtrations(max_dim+1);
        verbose = false;
        alpha_OFF(filename, max_dim, boundaries_triples, filtrations);
        std::vector<SparseMatrixInt> reindexed_boundaries(boundaries_triples.size());
        reindex_boundaries(boundaries_triples, reindexed_boundaries);
        set_boundaries_filtrations(reindexed_boundaries, filtrations);

    }

    Alpha::Alpha(std::vector<std::tuple<double,double,double>> points, int max_dim) : petls::Complex() { // call default constructor of parent class
        std::vector<std::vector<std::tuple<int,int,int>>> boundaries_triples(max_dim);
        std::vector<std::vector<filtration_type>> filtrations(max_dim+1);
        verbose = false;
        std::vector<Point> points_gudhi(points.size());
        std::cout << "number of points:" << points.size() << std::endl;
        for (int i = 0; i < points.size(); i++){
            points_gudhi[i] = Point(std::get<0>(points[i]),
                                    std::get<1>(points[i]),
                                    std::get<2>(points[i])
            );
        }
        alpha_points(points_gudhi, max_dim, boundaries_triples, filtrations);
        // alpha_OFF(filename, max_dim, boundaries_triples, filtrations);
        std::vector<SparseMatrixInt> reindexed_boundaries(boundaries_triples.size());
        reindex_boundaries(boundaries_triples, reindexed_boundaries);
        set_boundaries_filtrations(reindexed_boundaries, filtrations);
    }

    // functions for alpha complex
    Vertex_list from(const Cell_handle &ch)
    {
        Vertex_list the_list;
        for (auto i = 0; i < 4; i++)
        {
            the_list.push_back(ch->vertex(i));
        }
        return the_list;
    }
    Vertex_list from(const Facet &fct)
    {
        Vertex_list the_list;
        for (auto i = 0; i < 4; i++)
        {
            if (fct.second != i)
            {
                the_list.push_back(fct.first->vertex(i));
            }
        }
        return the_list;
    }
    Vertex_list from(const Edge &edg)
    {
        Vertex_list the_list;
        for (auto i = 0; i < 4; i++)
        {
            if ((edg.second == i) || (edg.third == i))
            {
                the_list.push_back(edg.first->vertex(i));
            }
        }
        return the_list;
    }
    Vertex_list from(const Alpha_shape_3::Vertex_handle &vh)
    {
        Vertex_list the_list;
        the_list.push_back(vh);
        return the_list;
    }


    void get_boundaries_and_filtrations(Simplex_tree simplex_tree, int dim_max,  std::vector<std::vector<std::tuple<int,int,int>>>& boundaries_triples,std::vector<std::vector<filtration_type>>& filtrations){
        // give each simplex a unique index (assign_key)
        // method from Gudhi rips_persistence_via_boundary_matrix.cpp
        int count = 0;
        for (auto simplex_handle : simplex_tree.filtration_simplex_range()){
            simplex_tree.assign_key(simplex_handle, count++);
        }

        // int max_dim = 3;

        //initialize list of boundaries
        for (int dim = 0; dim < dim_max; dim++){
            boundaries_triples[dim] = std::vector<std::tuple<int,int,int>>();
        }

        // initialize list of filtrations
        for (int dim = 0; dim <= dim_max; dim++){
            filtrations[dim] = std::vector<filtration_type>();
        }

        // iterate over all simplices to store their boundary information
        for (auto f_simplex : simplex_tree.filtration_simplex_range())
        {
            int dim_f = simplex_tree.dimension(f_simplex);

            // store filtration
            int f_key = simplex_tree.key(f_simplex);
            filtrations[dim_f].push_back(simplex_tree.filtration(f_simplex));

            // loop over the boundaries
            
            int sign = 1 - 2 * (dim_f % 2); // from Gudhi Persistent_cohomology.h
            for (auto b_simplex : simplex_tree.boundary_simplex_range(f_simplex))
            {

                boundaries_triples[dim_f-1].push_back(std::make_tuple(simplex_tree.key(b_simplex),
                                                                    f_key,
                                                                    sign));
                sign = -sign;
            }
        }
        
    }

    void alpha_points(std::vector<Point>& points, int dim_max, std::vector<std::vector<std::tuple<int,int,int>>>& boundaries_triples,std::vector<std::vector<filtration_type>>& filtrations){
        // alpha shape construction from points. CGAL has a strange behavior in REGULARIZED mode.
        Alpha_shape_3 as(points.begin(), points.end(), 0, Alpha_shape_3::GENERAL);

        // filtration with alpha values from alpha shape
        std::vector<Object> the_objects;
        std::vector<Alpha_value_type> the_alpha_values;
        Dispatch disp = CGAL::dispatch_output<Object, Alpha_value_type>(std::back_inserter(the_objects),
                                                                        std::back_inserter(the_alpha_values));
        as.filtration_with_alpha_values(disp);

        Alpha_shape_3::size_type count_vertices = 0;
        Alpha_shape_3::size_type count_edges = 0;
        Alpha_shape_3::size_type count_facets = 0;
        Alpha_shape_3::size_type count_cells = 0;
        // Loop on objects vector
        Vertex_list vertex_list;
        Simplex_tree simplex_tree;
        Alpha_shape_simplex_tree_map map_cgal_simplex_tree;
        std::vector<Alpha_value_type>::iterator the_alpha_value_iterator = the_alpha_values.begin();
        for (auto object_iterator : the_objects)
        {
            // Retrieve Alpha shape vertex list from object
            if (const Cell_handle *cell = CGAL::object_cast<Cell_handle>(&object_iterator))
            {
                vertex_list = from(*cell);
                count_cells++;
            }
            else if (const Facet *facet = CGAL::object_cast<Facet>(&object_iterator))
            {
                vertex_list = from(*facet);
                count_facets++;
            }
            else if (const Edge *edge = CGAL::object_cast<Edge>(&object_iterator))
            {
                vertex_list = from(*edge);
                count_edges++;
            }
            else if (const Alpha_shape_3::Vertex_handle *vertex =
                        CGAL::object_cast<Alpha_shape_3::Vertex_handle>(&object_iterator))
            {
                count_vertices++;
                vertex_list = from(*vertex);
            }
            // Construction of the vector of simplex_tree vertex from list of alpha_shapes vertex
            Simplex_tree_vector_vertex the_simplex_tree;
            for (auto the_alpha_shape_vertex : vertex_list)
            {
                Alpha_shape_simplex_tree_map::iterator the_map_iterator = map_cgal_simplex_tree.find(the_alpha_shape_vertex);
                if (the_map_iterator == map_cgal_simplex_tree.end())
                {
                    // alpha shape not found
                    Simplex_tree_vertex vertex = map_cgal_simplex_tree.size();

                    the_simplex_tree.push_back(vertex);
                    map_cgal_simplex_tree.insert(Alpha_shape_simplex_tree_pair(the_alpha_shape_vertex, vertex));
                }
                else
                {
                    // alpha shape found
                    Simplex_tree_vertex vertex = the_map_iterator->second;

                    the_simplex_tree.push_back(vertex);
                }
            }
            // Construction of the simplex_tree

            simplex_tree.insert_simplex(the_simplex_tree, std::sqrt(*the_alpha_value_iterator));
            if (the_alpha_value_iterator != the_alpha_values.end())
                ++the_alpha_value_iterator;
            else
                std::cerr << "This shall not happen" << std::endl;
        }
        get_boundaries_and_filtrations(simplex_tree, dim_max, boundaries_triples, filtrations);
        
    }
    

    void alpha_OFF(const char* filename, int dim_max, std::vector<std::vector<std::tuple<int,int,int>>>& boundaries_triples,std::vector<std::vector<filtration_type>>& filtrations){
        // Read points from file
        std::string offInputFile(filename);
        // Read the OFF file (input file name given as parameter) and triangulate points
        Gudhi::Points_3D_off_reader<Point> off_reader(offInputFile);
        // Check the read operation was correct
        if (!off_reader.is_valid())
        {
            std::cerr << "Unable to read file " << filename << std::endl;
            return;
        }
        // Retrieve the triangulation
        std::vector<Point> points = off_reader.get_point_cloud();
        alpha_points(points, dim_max, boundaries_triples, filtrations);
    }
}



#else // PETLS_USE_ALPHA_COMPLEX not defined
namespace petls{
    Alpha::Alpha(const char* filename, int max_dim){
        std::cout << "PETLS_USE_ALPHA_COMPLEX not defined; Alpha complex not being used." << std::endl;
    } 
    Alpha::Alpha(std::vector<std::tuple<double,double,double>> points, int max_dim){
        std::cout << "PETLS_USE_ALPHA_COMPLEX not defined; Alpha complex not being used." << std::endl;
    } 
}
#endif // PETLS_USE_ALPHA_COMPLEX not defined