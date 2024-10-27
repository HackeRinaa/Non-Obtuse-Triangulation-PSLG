#include "../HeaderFiles/triangulation_utils.hpp"
#include <iostream>
#include <cmath>
#include <vector>
#include <nlohmann/json.hpp>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/squared_distance_2.h>
#include <algorithm> // for std::clamp

// Function to parse the JSON input
void parseInput(const nlohmann::json &input, std::vector<Point> &points, std::vector<std::pair<int, int>> &constraints)
{
    if (input.contains("points_x") && input.contains("points_y"))
    {
        const auto &points_x = input["points_x"];
        const auto &points_y = input["points_y"];

        if (points_x.size() != input["num_points"] || points_y.size() != input["num_points"])
        {
            throw std::runtime_error("Number of points does not match the specified count.");
        }

        for (size_t i = 0; i < points_x.size(); ++i)
        {
            points.emplace_back(points_x[i].get<double>(), points_y[i].get<double>());
        }
    }
    else
    {
        throw std::runtime_error("Missing 'points_x' or 'points_y' in JSON.");
    }

    if (input.contains("additional_constraints"))
    {
        const auto &additional_constraints = input["additional_constraints"];
        if (!additional_constraints.empty())
        {
            for (const auto &constraint : additional_constraints)
            {
                if (constraint.size() == 2)
                {
                    int index1 = constraint[0].get<int>();
                    int index2 = constraint[1].get<int>();

                    if (index1 < points.size() && index2 < points.size())
                    {
                        constraints.emplace_back(index1, index2);
                    }
                    else
                    {
                        throw std::out_of_range("Constraint indices are out of bounds.");
                    }
                }
                else
                {
                    throw std::runtime_error("Invalid constraint format; each constraint must have exactly two indices.");
                }
            }
        }
        else
        {
            std::cout << "No additional constraints provided." << std::endl;
        }
    }
}

double approximate_angle_2D(const Point &p1, const Point &p2, const Point &p3)
{
    double a2 = CGAL::squared_distance(p2, p3);
    double b2 = CGAL::squared_distance(p1, p3);
    double c2 = CGAL::squared_distance(p1, p2);

    double b = std::sqrt(b2);
    double c = std::sqrt(c2);

    double cosine_angle = (b2 + c2 - a2) / (2 * b * c);
    cosine_angle = boost::algorithm::clamp(cosine_angle, -1.0, 1.0);

    double angle_radians = std::acos(cosine_angle);
    return angle_radians * (180.0 / M_PI);
}

bool isObtuseAngle(const Point &a, const Point &b, const Point &c)
{
    double ab2 = CGAL::squared_distance(a, b);
    double bc2 = CGAL::squared_distance(b, c);
    double ca2 = CGAL::squared_distance(c, a);

    return (ab2 + bc2 < ca2) || (ab2 + ca2 < bc2) || (bc2 + ca2 < ab2);
}

bool isObtuse(CDT::Face_handle face)
{
    if (face == nullptr)
        return false;
    return isObtuseAngle(face->vertex(0)->point(), face->vertex(1)->point(), face->vertex(2)->point()) ||
           isObtuseAngle(face->vertex(1)->point(), face->vertex(2)->point(), face->vertex(0)->point()) ||
           isObtuseAngle(face->vertex(2)->point(), face->vertex(0)->point(), face->vertex(1)->point());
}

Point obtuseAngleProjection(const Point &obtuse_vertex, const Point &p1, const Point &p2)
{
    Vector vec = p2 - p1;
    Vector dir = obtuse_vertex - p1;
    FT t = (dir * vec) / (vec * vec);
    Point projection = p1 + t * vec;
    return projection;
}

Point insertSteinerAtMedian(const Point &p1, const Point &p2)
{
    return CGAL::midpoint(p1, p2);
}

Point insertSteinerAtCircumcenter(const Point &p1, const Point &p2, const Point &p3)
{
    return CGAL::circumcenter(p1, p2, p3);
}

Point insertSteinerForPolygon(const std::vector<Point> &polygon_points)
{
    if (polygon_points.empty())
    {
        throw std::runtime_error("Polygon points cannot be empty.");
    }
    FT x_sum = 0, y_sum = 0;
    for (const auto &point : polygon_points)
    {
        x_sum += point.x();
        y_sum += point.y();
    }
    FT n = static_cast<FT>(polygon_points.size());
    return Point(x_sum / n, y_sum / n);
}

void insertSteinerPoint(CDT &cdt, CDT::Face_handle face)
{
    if (face == nullptr)
    {
        std::cerr << "Error: Invalid face provided for Steiner point insertion." << std::endl;
        return;
    }

    Point p1 = face->vertex(0)->point();
    Point p2 = face->vertex(1)->point();
    Point p3 = face->vertex(2)->point();

    int obtuse_index = -1;
    if (isObtuseAngle(p1, p2, p3))
        obtuse_index = 0;
    else if (isObtuseAngle(p2, p3, p1))
        obtuse_index = 1;
    else if (isObtuseAngle(p3, p1, p2))
        obtuse_index = 2;

    Point steiner_point;
    if (obtuse_index != -1)
    {
        switch (obtuse_index)
        {
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
    }
    else
    {
        steiner_point = insertSteinerAtCircumcenter(p1, p2, p3);
    }

    CDT::Vertex_handle vh = cdt.insert(steiner_point, face);
    if (vh == nullptr)
    {
        std::cerr << "Error: Steiner point insertion failed." << std::endl;
    }
    else
    {
        std::cout << "Inserted Steiner point at: " << vh->point() << std::endl;
    }
}

void applyPolygonalTreatment(CDT &cdt, CDT::Face_handle face1, CDT::Face_handle face2)
{
    std::vector<Point> polygon_points;
    if (face1 == nullptr || face2 == nullptr)
    {
        std::cerr << "Error: One of the faces is null or invalid." << std::endl;
        return;
    }
    for (int i = 0; i < 3; ++i)
    {
        polygon_points.push_back(face1->vertex(i)->point());
        polygon_points.push_back(face2->vertex(i)->point());
    }

    Point steiner_point = insertSteinerForPolygon(polygon_points);

    CDT::Vertex_handle shared_v1 = nullptr;
    CDT::Vertex_handle shared_v2 = nullptr;
    int shared_edge_index = -1;

    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            if (face1->vertex(i) == face2->vertex(j))
            {
                if (shared_v1 == nullptr)
                {
                    shared_v1 = face1->vertex(i);
                }
                else if (shared_v2 == nullptr)
                {
                    shared_v2 = face1->vertex(i);
                    shared_edge_index = i;
                }
            }
        }
    }

    if (shared_v1 == nullptr || shared_v2 == nullptr || shared_edge_index == -1)
    {
        std::cerr << "Error: No shared edge found between the faces." << std::endl;
        return;
    }

    if (cdt.is_constrained(CDT::Edge(face1, shared_edge_index)))
    {
        cdt.remove_constraint(face1, shared_edge_index);
    }
    else
    {
        std::cerr << "Error: The shared edge is not constrained." << std::endl;
    }
}

// Function to check if an edge flip is beneficial
bool isEdgeFlipBeneficial(CDT &cdt, CDT::Face_handle &face, int edge_index)
{
    CDT::Face_handle neighbor = face->neighbor(edge_index);

    if (!cdt.is_infinite(neighbor))
    {
        // Get the vertices of the edge
        CDT::Vertex_handle v1 = face->vertex((edge_index + 1) % 3);
        CDT::Vertex_handle v2 = face->vertex((edge_index + 2) % 3);

        // Get the opposite vertices in both faces
        CDT::Vertex_handle opp1 = face->vertex(edge_index);
        CDT::Vertex_handle opp2 = neighbor->vertex(neighbor->index(face));

        // Calculate the angles in the new triangles formed after a flip
        K::FT angle1 = approximate_angle_2D(v1->point(), opp2->point(), v2->point());
        K::FT angle2 = approximate_angle_2D(v2->point(), opp1->point(), v1->point());

        if (angle1 < 90 && angle2 < 90)
        {
            return true;
        }
    }
    return false;
}

bool applyEdgeFlip(CDT &cdt, CDT::Face_handle &face)
{
    // Try flipping each edge of the triangle and see if it improves the triangulation
    for (int i = 0; i < 3; ++i)
    {
        CDT::Face_handle neighbor = face->neighbor(i); // Get the adjacent face across edge `i`
        if (!cdt.is_infinite(neighbor))
        {
            // Flip the edge and check if the obtuse angle is resolved
            if (isEdgeFlipBeneficial(cdt, face, i))
            {
                cdt.flip(face, i);
                return true;
            }
        }
    }

    return false;
}

void writeOutput(const CDT &cdt, const std::vector<Point> &steinerPoints, nlohmann::json &output)
{
    std::vector<double> steiner_x, steiner_y;
    for (const auto &point : steinerPoints)
    {
        steiner_x.push_back(CGAL::to_double(point.x()));
        steiner_y.push_back(CGAL::to_double(point.y()));
    }
    output["steiner_points_x"] = steiner_x;
    output["steiner_points_y"] = steiner_y;

    std::vector<int> face_indices;
    for (auto it = cdt.finite_faces_begin(); it != cdt.finite_faces_end(); ++it)
    {
        for (int i = 0; i < 3; ++i)
        {
            face_indices.push_back(it->vertex(i)->info());
        }
    }
    output["face_indices"] = face_indices;
}
