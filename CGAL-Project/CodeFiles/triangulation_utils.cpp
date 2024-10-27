#include "../HeaderFiles/triangulation_utils.hpp"
#include <iostream>
#include <cmath>
#include <vector>
#include <CGAL/squared_distance_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>

// Function to parse the JSON input
void parseInput(const nlohmann::json& input, std::vector<Point>& points, std::vector<std::pair<int, int>>& constraints) {
    // Check if 'points_x' and 'points_y' keys exist
    if (input.contains("points_x") && input.contains("points_y")) {
        const auto& points_x = input["points_x"];
        const auto& points_y = input["points_y"];
        
        // Ensure that the number of points matches the specified number
        if (points_x.size() != input["num_points"] || points_y.size() != input["num_points"]) {
            throw std::runtime_error("Number of points does not match the specified count.");
        }

        // Read the points from the JSON
        for (size_t i = 0; i < points_x.size(); ++i) {
            points.emplace_back(points_x[i].get<double>(), points_y[i].get<double>());
        }
    } else {
        throw std::runtime_error("Missing 'points_x' or 'points_y' in JSON.");
    }

    // Check if 'additional_constraints' key exists
    if (input.contains("additional_constraints") ) {
        const auto& additional_constraints = input["additional_constraints"];

        // Read the constraints from the JSON
        for (const auto& constraint : additional_constraints) {
            if (constraint.size() == 2) {
                int index1 = constraint[0].get<int>();
                int index2 = constraint[1].get<int>();
                
                // Ensure the indices are valid
                if (index1 < points.size() && index2 < points.size()) {
                    constraints.emplace_back(index1, index2);
                } else {
                    throw std::out_of_range("Constraint indices are out of bounds.");
                }
            } else {
                throw std::runtime_error("Invalid constraint format; each constraint must have exactly two indices.");
            }
        }
    }
}

double approximate_angle_2D(const Point& p1, const Point& p2, const Point& p3) {
    // Calculate squared distances for the triangle sides
    double a2 = CGAL::squared_distance(p2, p3);
    double b2 = CGAL::squared_distance(p1, p3);
    double c2 = CGAL::squared_distance(p1, p2);
    
    // Calculate the lengths of the sides
    double b = std::sqrt(b2);
    double c = std::sqrt(c2);

    // Apply the law of cosines to find the angle at p2
    double cosine_angle = (b2 + c2 - a2) / (2 * b * c);
    cosine_angle =boost::algorithm::clamp(cosine_angle, -1.0, 1.0);  // Ensure the value stays within [-1, 1]

    // Calculate the angle in radians and convert it to degrees
    double angle_radians = std::acos(cosine_angle);
    double angle_degrees = angle_radians * (180.0 / M_PI);

    return angle_degrees;
}


// Helper function to check for obtuse angles
bool isObtuseAngle(const Point& a, const Point& b, const Point& c) {
    // Calculate the squared length of the sides
    double ab2 = CGAL::squared_distance(a, b);
    double bc2 = CGAL::squared_distance(b, c);
    double ca2 = CGAL::squared_distance(c, a);

    // Apply the triangle inequality to check if the angle at b is obtuse
    return (ab2 + bc2 < ca2) || (ab2 + ca2 < bc2) || (bc2 + ca2 < ab2);
}

// Check if the triangle has an obtuse angle
bool isObtuse(CDT::Face_handle face) {
    if (face == nullptr) return false;

    // Check if any of the three angles in the face is obtuse
    return isObtuseAngle(face->vertex(0)->point(), face->vertex(1)->point(), face->vertex(2)->point()) ||
           isObtuseAngle(face->vertex(1)->point(), face->vertex(2)->point(), face->vertex(0)->point()) ||
           isObtuseAngle(face->vertex(2)->point(), face->vertex(0)->point(), face->vertex(1)->point());
}



// 1. Obtuse Angle Projection Strategy
Point obtuseAngleProjection(const Point& obtuse_vertex, const Point& p1, const Point& p2) {
    // Project obtuse_vertex onto line segment (p1, p2)
    Vector vec = p2 - p1;
    Vector dir = obtuse_vertex - p1;
    FT t = (dir * vec) / (vec * vec);
    Point projection = p1 + t * vec;
    return projection;
}

// 2. Median Insertion Strategy
Point insertSteinerAtMedian(const Point& p1, const Point& p2) {
    // Insert at the midpoint of the longest edge
    return CGAL::midpoint(p1, p2);
}

// 3. Circumcenter Insertion Strategy
Point insertSteinerAtCircumcenter(const Point& p1, const Point& p2, const Point& p3) {
    // Insert at the circumcenter of the triangle
    Point circumcenter = CGAL::circumcenter(p1, p2, p3);
    return circumcenter;
}

// 4. Polygonal Treatment Strategy
Point insertSteinerForPolygon(const std::vector<Point>& polygon_points) {
    // Check if the input polugon points are empty
    if (polygon_points.empty()) {
        throw std::runtime_error("Polygon points cannot be empty.");
    }

    // Sum of x and y coordinates
    FT x_sum = 0, y_sum = 0;
    for (const auto& point : polygon_points) {
        x_sum += point.x();
        y_sum += point.y();
    }

    // Calculate the centroid of the polygon and return it as the Steiner point
    FT n = static_cast<FT>(polygon_points.size());
    return Point(x_sum / n, y_sum / n);
}

// Main Steiner insertion function
void insertSteinerPoint(CDT& cdt, CDT::Face_handle face) {
    if (face == nullptr) {
        std::cerr << "Error: Invalid face provided for Steiner point insertion." << std::endl;
        return;
    }

    // Extract the three vertices of the face
    Point p1 = face->vertex(0)->point();
    Point p2 = face->vertex(1)->point();
    Point p3 = face->vertex(2)->point();

    // Check if the face contains an obtuse angle and determine the obtuse vertex
    int obtuse_index = -1;
    if (isObtuseAngle(p1, p2, p3)) obtuse_index = 0;
    else if (isObtuseAngle(p2, p3, p1)) obtuse_index = 1;
    else if (isObtuseAngle(p3, p1, p2)) obtuse_index = 2;

    Point steiner_point;
    if (obtuse_index != -1) {
        // There is an obtuse angle in the face; choose the obtuse angle projection strategy
        switch (obtuse_index) {
            case 0:
                steiner_point = obtuseAngleProjection(p1, p2, p3);
                break;
            case 1:
                steiner_point = obtuseAngleProjection(p2, p3, p1);
                break;
            case 2:
                steiner_point = obtuseAngleProjection(p3, p1, p2);
                break;
        }
    } else {
        // If no obtuse angle is found, fall back to a different insertion strategy, e.g., circumcenter
        steiner_point = insertSteinerAtCircumcenter(p1, p2, p3);
    }

    // Insert the Steiner point into the triangulation
    CDT::Vertex_handle vh = cdt.insert(steiner_point, face);
    if (vh == nullptr) {
        std::cerr << "Error: Steiner point insertion failed." << std::endl;
    } else {
        std::cout << "Inserted Steiner point at: " << vh->point() << std::endl;
    }
}


void applyPolygonalTreatment(CDT& cdt, CDT::Face_handle face1, CDT::Face_handle face2) {
    std::vector<Point> polygon_points;

    if (face1 == nullptr || face2 == nullptr) {
        std::cerr << "Error: One of the faces is null or invalid." << std::endl;
        return;
    }

    // Collect all unique vertices from both faces
    for (int i = 0; i < 3; ++i) {
        polygon_points.push_back(face1->vertex(i)->point());
        polygon_points.push_back(face2->vertex(i)->point());
    }

    Point steiner_point = insertSteinerForPolygon(polygon_points);

    CDT::Vertex_handle shared_v1 = nullptr;
    CDT::Vertex_handle shared_v2 = nullptr;
    int shared_edge_index = -1;

    // Identify the shared edge between face1 and face2
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            if (face1->vertex(i) == face2->vertex(j)) {
                if (shared_v1 == nullptr) {
                    shared_v1 = face1->vertex(i);
                } else if (shared_v2 == nullptr) {
                    shared_v2 = face1->vertex(i);
                    shared_edge_index = i;
                }
            }
        }
    }

    if (shared_v1 == nullptr || shared_v2 == nullptr || shared_edge_index == -1) {
        std::cerr << "Error: No shared edge found between the faces." << std::endl;
        return;
    }

    // Check if the shared edge is constrained and remove it if so
    if (cdt.is_constrained(CDT::Edge(face1, shared_edge_index))) {
        cdt.remove_constraint(face1, shared_edge_index);
    } else {
        std::cerr << "Error: The shared edge is not constrained." << std::endl;
        return;
    }

    cdt.insert(steiner_point);
}


// Function to get the adjacent face of a given face across a specific edge
Face_handle getAdjacentFace(Face_handle face, int edge_index) {
    // `face->neighbor(edge_index)` gets the adjacent face across the `edge_index` edge
    return face->neighbor(edge_index);
}

// Function to find an adjacent face that is obtuse
CDT::Face_handle findObtuseAdjacentFace(CDT& cdt, CDT::Face_handle face) {
    for (int i = 0; i < 3; ++i) {
        CDT::Face_handle neighbor = face->neighbor(i);

        if (!cdt.is_infinite(neighbor)) {
            if (isObtuse(neighbor)) {
                return neighbor;
            }
        }
    }
    return nullptr;  // Return nullptr if no obtuse adjacent face is found
}



// Function to check if an edge flip is beneficial
bool isEdgeFlipBeneficial(CDT& cdt, CDT::Face_handle& face, int edge_index) {
    // Get the neighboring face across the edge
    CDT::Face_handle neighbor = face->neighbor(edge_index);

    // Ensure the neighbor is not infinite
    if (!cdt.is_infinite(neighbor)) {
        // Get the vertices of the edge
        CDT::Vertex_handle v1 = face->vertex((edge_index + 1) % 3);
        CDT::Vertex_handle v2 = face->vertex((edge_index + 2) % 3);

        // Get the opposite vertices in both faces
        CDT::Vertex_handle opp1 = face->vertex(edge_index);
        CDT::Vertex_handle opp2 = neighbor->vertex(neighbor->index(face));

        // Calculate the angles in the new triangles formed after a flip
        K::FT angle1 = approximate_angle_2D(v1->point(), opp2->point(), v2->point());
        K::FT angle2 = approximate_angle_2D(v2->point(), opp1->point(), v1->point());

        // Check if the new angles are less than 90 degrees
        if (angle1 < 90 && angle2 < 90) {
            return true;  // Edge flip is beneficial
        }
    }
    return false;
}


bool applyEdgeFlip(CDT& cdt, CDT::Face_handle& face) {
    // Try flipping each edge of the triangle and see if it improves the triangulation
    for (int i = 0; i < 3; ++i) {
        CDT::Face_handle neighbor = face->neighbor(i);  // Get the adjacent face across edge `i`
        if (!cdt.is_infinite(neighbor)) {
            // Flip the edge and check if the obtuse angle is resolved
            if (isEdgeFlipBeneficial(cdt, face, i)) {
                cdt.flip(face, i);  // Perform edge flip
                return true;  // Flip was successful and reduced the obtuse angle
            }
        }
    }

    return false;  // No beneficial flip was found
}


// Function to write output to a file
void writeOutput(const CDT& cdt, const std::vector<Point>& steinerPoints, nlohmann::json& output) {
    output["faces"] = nlohmann::json::array();
    output["steinerPoints"] = nlohmann::json::array();

    // Write the faces
    for (auto fit = cdt.finite_faces_begin(); fit != cdt.finite_faces_end(); ++fit) {
        std::array<double, 2> p1 = {fit->vertex(0)->point().x(), fit->vertex(0)->point().y()};
        std::array<double, 2> p2 = {fit->vertex(1)->point().x(), fit->vertex(1)->point().y()};
        std::array<double, 2> p3 = {fit->vertex(2)->point().x(), fit->vertex(2)->point().y()};
        output["faces"].push_back({p1, p2, p3});
    }

    // Write the Steiner points
    for (const auto& point : steinerPoints) {
        output["steinerPoints"].push_back({point.x(), point.y()});
    }
}

