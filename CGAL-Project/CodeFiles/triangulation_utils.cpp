#include "../HeaderFiles/triangulation_utils.h"
#include <cmath>
#include <iostream>

// Check if a triangle is obtuse
bool isObtuse(const CDT::Face_handle &face)
{
    Point p0 = face->vertex(0)->point();
    Point p1 = face->vertex(1)->point();
    Point p2 = face->vertex(2)->point();

    double a2 = CGAL::squared_distance(p1, p2);
    double b2 = CGAL::squared_distance(p0, p2);
    double c2 = CGAL::squared_distance(p0, p1);

    if ((a2 + b2) < c2)
        return true; // Angle at p0 is obtuse
    if ((a2 + c2) < b2)
        return true; // Angle at p1 is obtuse
    if ((b2 + c2) < a2)
        return true; // Angle at p2 is obtuse

    return false; // No obtuse angles
}

// Calculate the Steiner point for an obtuse triangle
Point calculateSteinerPoint(const CDT::Face_handle &face)
{
    Point p0 = face->vertex(0)->point();
    Point p1 = face->vertex(1)->point();
    Point p2 = face->vertex(2)->point();

    double a2 = CGAL::squared_distance(p1, p2);
    double b2 = CGAL::squared_distance(p0, p2);
    double c2 = CGAL::squared_distance(p0, p1);

    // Find the longest edge and return its midpoint
    if (a2 >= b2 && a2 >= c2)
    {
        return CGAL::midpoint(p1, p2);
    }
    else if (b2 >= a2 && b2 >= c2)
    {
        return CGAL::midpoint(p0, p2);
    }
    else
    {
        return CGAL::midpoint(p0, p1);
    }
}

// Attempt to flip diagonals to reduce obtuse angles
bool tryDiagonalFlip(CDT &cdt, CDT::Face_handle face)
{
    for (int i = 0; i < 3; ++i)
    {
        if (cdt.is_flipable(face, i))
        {
            cdt.flip(face, i); // Correct usage of flip method
            std::cout << "Diagonal flip applied." << std::endl;
            return true; // A flip was successful
        }
    }
    return false; // No flip was possible
}

void insertSteinerPoints(CDT &cdt)
{
    bool has_obtuse_angle = true;
    int max_iterations = 1000; // Limit to avoid infinite loops
    int iteration = 0;
    int obtuse_angles = 0;             // Count of obtuse angles found in each iteration
    static int prev_obtuse_angles = 0; // Declare prev_obtuse_angles as a static variable

    while (has_obtuse_angle && iteration < max_iterations)
    {
        has_obtuse_angle = false;
        obtuse_angles = 0; // Reset count for each iteration

        std::vector<Point> steinerPoints;

        // Iterate over all finite faces in the triangulation
        for (auto face = cdt.finite_faces_begin(); face != cdt.finite_faces_end(); ++face)
        {
            if (!tryDiagonalFlip(cdt, face))
            { // Correct call without reference
                if (isObtuse(face))
                {
                    Point steinerPoint = calculateSteinerPoint(face);
                    steinerPoints.push_back(steinerPoint); // Collect Steiner points
                    has_obtuse_angle = true;               // Flag to keep iterating
                    obtuse_angles++;                       // Increment count of obtuse angles
                    std::cout << "Found obtuse angle, inserting Steiner point: " << steinerPoint << std::endl;
                }
            }
        }

        // Insert Steiner points in a separate loop
        for (const auto &point : steinerPoints)
        {
            cdt.insert(point);
        }

        iteration++;
        std::cout << "Iteration " << iteration << ": " << steinerPoints.size() << " Steiner points inserted. Found " << obtuse_angles << " obtuse angles.\n";

        // Check if the number of obtuse angles has decreased
        if (obtuse_angles == 0 || obtuse_angles == prev_obtuse_angles)
        {
            break; // Stop iterating if no obtuse angles are found or the count hasn't decreased
        }
        prev_obtuse_angles = obtuse_angles; // Store the previous count for comparison
    }

    if (iteration >= max_iterations)
    {
        std::cerr << "Warning: Max iterations reached. The triangulation may still contain obtuse angles." << std::endl;
    }
}
