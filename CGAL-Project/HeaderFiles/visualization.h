#ifndef VISUALIZATION_H
#define VISUALIZATION_H

#include <QMainWindow>
#include <QTextEdit>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Constrained_Delaunay_triangulation_2<K> CDT;
typedef CDT::Point Point;

class VisualizationWindow : public QMainWindow
{
    Q_OBJECT

public:
    VisualizationWindow(QWidget *parent = nullptr);
    void displayTriangulation(const CDT &cdt);

private:
    QTextEdit *textEdit;
};

#endif // VISUALIZATION_H
