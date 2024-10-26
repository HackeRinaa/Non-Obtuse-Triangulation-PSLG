#include "../HeaderFiles/triangulation_utils.h"
#include <cmath>
#include <iostream>
#include <set>

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
            cdt.flip(face, i); // Apply flip
            std::cout << "Diagonal flip applied." << std::endl;
            return true; // A flip was successful
        }
    }
    return false; // No flip was possible
}

// Attempt to flip diagonals more aggressively to reduce obtuse angles
bool tryEnhancedDiagonalFlip(CDT &cdt, CDT::Face_handle face)
{
    bool flip_successful = false;
    // Try flipping each edge of the triangle
    for (int i = 0; i < 3; ++i)
    {
        if (cdt.is_flipable(face, i))
        {
            cdt.flip(face, i);
            std::cout << "Enhanced diagonal flip applied.\n";
            flip_successful = true;
        }
    }
    return flip_successful;
}

// Updated insertSteinerPoints function with enhanced strategies
// Try flipping edges on neighboring faces of the obtuse triangle
bool tryExtendedFlip(CDT &cdt, CDT::Face_handle face)
{
    bool flip_successful = false;
    for (int i = 0; i < 3; ++i)
    {
        // Check if edge is flippable and perform the flip if it might help
        if (cdt.is_flipable(face, i))
        {
            cdt.flip(face, i);
            std::cout << "Extended flip applied on face.\n";
            flip_successful = true;
        }

        // Attempt flipping on adjacent faces as well
        auto neighbor = face->neighbor(i);
        for (int j = 0; j < 3; ++j)
        {
            if (cdt.is_flipable(neighbor, j))
            {
                cdt.flip(neighbor, j);
                std::cout << "Flip applied on adjacent face.\n";
                flip_successful = true;
            }
        }
    }
    return flip_successful;
}

// Attempt alternate Steiner point positions to minimize insertions
Point findAlternateSteinerPoint(const CDT::Face_handle &face)
{
    Point p0 = face->vertex(0)->point();
    Point p1 = face->vertex(1)->point();
    Point p2 = face->vertex(2)->point();

    // Calculate squared edge lengths
    double a2 = CGAL::squared_distance(p1, p2);
    double b2 = CGAL::squared_distance(p0, p2);
    double c2 = CGAL::squared_distance(p0, p1);

    // Try midpoints of different edges to find the best Steiner point
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

// Updated insertSteinerPoints function with extended flipping and alternate Steiner point positions
void insertSteinerPoints(CDT &cdt)
{
    bool changes_made = true;
    int max_iterations = 1000;
    int iteration = 0;
    std::set<CDT::Face_handle> processed_faces;

    while (changes_made && iteration < max_iterations)
    {
        changes_made = false;
        std::vector<Point> steinerPoints;

        for (auto face = cdt.finite_faces_begin(); face != cdt.finite_faces_end(); ++face)
        {
            if (processed_faces.count(face) == 0 && isObtuse(face))
            {
                // Try extended flip strategy on this face and its neighbors
                bool flip_successful = tryExtendedFlip(cdt, face);

                if (!flip_successful)
                {
                    // If no flip resolved the obtuse angle, consider an alternate Steiner point
                    Point steinerPoint = findAlternateSteinerPoint(face);
                    if (std::find(steinerPoints.begin(), steinerPoints.end(), steinerPoint) == steinerPoints.end())
                    {
                        steinerPoints.push_back(steinerPoint);
                        std::cout << "Inserting alternate Steiner point: " << steinerPoint << "\n";
                    }
                    changes_made = true;
                }

                processed_faces.insert(face); // Mark face as processed
            }
        }

        // Insert Steiner points
        for (const auto &point : steinerPoints)
        {
            cdt.insert(point);
        }

        iteration++;
        std::cout << "Iteration " << iteration << ": " << steinerPoints.size() << " Steiner points inserted.\n";

        if (!changes_made && steinerPoints.empty())
        {
            break;
        }
    }

    if (iteration >= max_iterations)
    {
        std::cerr << "Warning: Max iterations reached. The triangulation may still contain obtuse angles.\n";
    }
}

// void insertSteinerPoints(CDT &cdt)
// {
//     bool changes_made = true;
//     int max_iterations = 1000; // Limit to avoid infinite loops
//     int iteration = 0;
//     std::set<CDT::Face_handle> processed_faces; // Track processed faces

//     while (changes_made && iteration < max_iterations)
//     {
//         changes_made = false;             // Reset changes for this iteration
//         std::vector<Point> steinerPoints; // List of Steiner points to add

//         // Iterate over all finite faces in the triangulation
//         for (auto face = cdt.finite_faces_begin(); face != cdt.finite_faces_end(); ++face)
//         {
//             if (processed_faces.count(face) == 0) // Skip faces already processed
//             {
//                 if (isObtuse(face))
//                 {
//                     bool flip_successful = tryDiagonalFlip(cdt, face); // Try flipping the diagonal
//                     if (!flip_successful)
//                     {
//                         // If no flip was possible, add a Steiner point
//                         Point steinerPoint = calculateSteinerPoint(face);
//                         steinerPoints.push_back(steinerPoint); // Collect Steiner points
//                         changes_made = true;                   // Indicate a change has been made
//                         std::cout << "Found obtuse angle, inserting Steiner point: " << steinerPoint << std::endl;
//                     }

//                     processed_faces.insert(face); // Mark face as processed
//                 }
//             }
//         }

//         // Insert Steiner points in a separate loop
//         for (const auto &point : steinerPoints)
//         {
//             cdt.insert(point);
//         }

//         iteration++;
//         std::cout << "Iteration " << iteration << ": " << steinerPoints.size() << " Steiner points inserted.\n";

//         // Stop if no changes were made in this iteration
//         if (!changes_made && steinerPoints.empty())
//         {
//             break;
//         }
//     }

//     if (iteration >= max_iterations)
//     {
//         std::cerr << "Warning: Max iterations reached. The triangulation may still contain obtuse angles." << std::endl;
//     }
// }