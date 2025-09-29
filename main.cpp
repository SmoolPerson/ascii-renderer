#include <iostream>
#include <ctime>
#include <cmath>
#include <thread>
#include <vector>
#include <cstdlib>
#include <algorithm>
#include <fstream>
double transl_distance;
// very simple container with easy string
class Point {
public:
    double y;
    double x;

    Point() {
        x = 0;
        y = 0;
    }

    Point(double inputx, double inputy) {
        x = inputx;
        y = inputy;
    }

    void print() {
        std::cout << "(" << x << ", " << y << ")" << std::endl;
    }

    Point scale(double scalar) {
        return Point(scalar * x, scalar * y);
    }

    Point add(Point p) {
        return Point(x + p.x, y + p.y);
    }

    Point negate() {
        return Point(-x, -y);
    }

    Point normalize() {
        double length = sqrt(x*x + y*y);
        if (length == 0) {
            return Point(0, 0);
        }
        return Point(x / length, y / length);
    }
};

class Triangle {
    public:
        Point A;
        Point B;
        Point C;
        std::string material = "q";
        Triangle(Point A, Point B, Point C) {
            this->A = A;
            this->B = B;
            this->C = C;
        }
        Triangle() {
            A = Point(0, 0);
            B = Point(0, 0);
            C = Point(0, 0);
        }

        bool point_in_triangle(Point p) {
            // https://blackpawn.com/texts/pointinpoly/
            // https://stackoverflow.com/questions/2049582/how-to-determine-if-a-point-is-in-a-2d-triangle
            Point p0 = A;
            Point p1 = B;
            Point p2 = C;
            // Explanation: this is the solution to this vector equationp = p0 + (p1 - p0) * s + (p2 - p0) * t
            // If you're confused, look at the blackpawn link as that is where this all came from
            // ChatGPT generated the solution so thanks OPenAI
            double s = ((p.x - p0.x)*(p2.y - p0.y) - (p2.x - p0.x)*(p.y - p0.y)) / 
           ((p1.x - p0.x)*(p2.y - p0.y) - (p2.x - p0.x)*(p1.y - p0.y));

            double t = ((p1.x - p0.x)*(p.y - p0.y) - (p.x - p0.x)*(p1.y - p0.y)) / 
           ((p1.x - p0.x)*(p2.y - p0.y) - (p2.x - p0.x)*(p1.y - p0.y));

           if (s < 0 || t < 0) {
            return false;
           }
           else if (s > 1 || t > 1) {
            return false;
           }
           else if (s + t > 1) {
            return false;
           }
           else {
            return true;
           }
        }

        bool point_on_triangle(Point p) {
            Point p0 = A;
            Point p1 = B;
            Point p2 = C;
            double s = ((p.x - p0.x)*(p2.y - p0.y) - (p2.x - p0.x)*(p.y - p0.y)) / 
           ((p1.x - p0.x)*(p2.y - p0.y) - (p2.x - p0.x)*(p1.y - p0.y));

            double t = ((p1.x - p0.x)*(p.y - p0.y) - (p.x - p0.x)*(p1.y - p0.y)) / 
           ((p1.x - p0.x)*(p2.y - p0.y) - (p2.x - p0.x)*(p1.y - p0.y));
            if (fabs(s) < 0.05 || fabs(t) < 0.05 || fabs(s + t - 1) < 0.05) {
                return true;
            }
            return false;
        }

        void print() {
            std::cout << "Printing triangle: " << std::endl;
            std::cout << "p1: ";
            A.print();
            std::cout << "p2: ";
            B.print();
            std::cout << "p3: ";
            C.print();
            std::cout << "\n" << std::endl;
        }
};

struct BoundingBox {
    double min_x, max_x;
    double min_y, max_y;
    Triangle triangle;
};

class Point3D {
    public:
        double x;
        double y;
        double z;
        Point3D() {
            x = 0;
            y = 0;
            z = 0;
        }
        Point3D(double newx, double newy, double newz) {
            x = newx;
            y = newy;
            z = newz;
        }
        void print() {
        std::cout << "(" << x << ", " << y << ", " << z << ")" << std::endl;
        }

        Point3D scale(double scalar) {
            return Point3D(scalar * x, scalar * y, scalar * z);
        }

        Point3D add(Point3D p) {
            return Point3D(x + p.x, y + p.y, z + p.z);
        }

        Point3D negate() {
            return Point3D(-x, -y, -z);
        }

        Point3D normalize() {
            double length = sqrt(x*x + y*y + z*z);
            if (length == 0) {
                return Point3D(1,0,0);
            }
            return Point3D(x / length, y / length, z / length);
        }

        Point3D cross(Point3D other) {
            return Point3D(y*other.z - z*other.y, z*other.x - x*other.z, x*other.y - y*other.x);
        }

        double dot(Point3D other) {
            return (x * other.x + y * other.y + z * other.z);
        }
};

class Triangle3D {
    public:
        Point3D A;
        Point3D B;
        Point3D C;
        std::string material = "";

        Triangle3D(Point3D A, Point3D B, Point3D C) {
            this->A = A;
            this->B = B;
            this->C = C;
        }
        Triangle3D() {
            A = Point3D(0, 0, 0);
            B = Point3D(0, 0, 0);
            C = Point3D(0, 0, 0);
        }

        void print() {
            std::cout << "Printing triangle: " << std::endl;
            std::cout << "p1: ";
            A.print();
            std::cout << "p2: ";
            B.print();
            std::cout << "p3: ";
            C.print();
            std::cout << "\n" << std::endl;
        }
};

// very simple matrix class just supporting matrix-vector multplication
class Matrix {
    public:
        double m_values[3][3];
        Matrix(double r1c1, double r1c2, double r1c3, double r2c1, double r2c2, double r2c3, double r3c1, double r3c2, double r3c3) {
            m_values[0][0] = r1c1;
            m_values[0][1] = r1c2;
            m_values[0][2] = r1c3;

            m_values[1][0] = r2c1;
            m_values[1][1] = r2c2;
            m_values[1][2] = r2c3;

            m_values[2][0] = r3c1;
            m_values[2][1] = r3c2;
            m_values[2][2] = r3c3;
        }

        Point3D multiply_with_vector(Point3D vector) {
            return Point3D(m_values[0][0] * vector.x + m_values[0][1] * 
                vector.y + m_values[0][2] * vector.z, 
                m_values[1][0] * vector.x + m_values[1][1] * 
                vector.y + m_values[1][2] * vector.z,
                m_values[2][0] * vector.x + m_values[2][1] * 
                vector.y + m_values[2][2] * vector.z);
        }
};

// a board for where points can be printed
class Board {
    public:
        BoundingBox tri_bounding_box;

        Board(double min_x, double min_y, double max_x, double max_y) {
            tri_bounding_box.min_x = min_x;
            tri_bounding_box.max_x = max_x;
            tri_bounding_box.max_y = max_y;
            tri_bounding_box.min_y = min_y;
        }

        void print_board(std::vector<Triangle> triangle, int length) {
            std::cout << "\033[H\033[2J"; // clear the console

            int width = tri_bounding_box.max_x - tri_bounding_box.min_x + 1;
            int height = tri_bounding_box.max_y - tri_bounding_box.min_y + 1; // index at an offset to prevent negative indices
            std::vector<std::vector<std::string>> framebuffer(width, std::vector<std::string>(height, ".."));
            
            for (size_t l = 0; l < triangle.size(); l++) {
                BoundingBox temp;

                // compute bounding box then loop over that
                temp.max_x = std::max(triangle[l].A.x, std::max(triangle[l].B.x, triangle[l].C.x)) - tri_bounding_box.min_x;
                temp.max_y = std::max(triangle[l].A.y, std::max(triangle[l].B.y, triangle[l].C.y)) - tri_bounding_box.min_y;
                temp.min_x = std::min(triangle[l].A.x, std::min(triangle[l].B.x, triangle[l].C.x)) - tri_bounding_box.min_x;
                temp.min_y = std::min(triangle[l].A.y, std::min(triangle[l].B.y, triangle[l].C.y)) - tri_bounding_box.min_y;;
                
                // clamp to make sure it isn't out of the screen
                int minX = std::max(0, (int) temp.min_x);
                int maxX = std::min(width - 1, (int) temp.max_x);
                int minY = std::max(0, (int) temp.min_y);
                int maxY = std::min(height - 1, (int) temp.max_y);
                for (int x = minX; x <= maxX; x++) {
                    for (int y = minY; y <= maxY; y++) {
                        // undo offset for point in triangle test
                        if (triangle[l].point_in_triangle(Point(x + tri_bounding_box.min_x, y + tri_bounding_box.min_y))) {
                            framebuffer[x][y] = triangle[l].material + triangle[l].material;
                        }
                    }
                }
            }
        
            // Iterate down since in the console things are printed from up to down
            for (int y = tri_bounding_box.max_y; y >= tri_bounding_box.min_y; --y) {
                for (int x = tri_bounding_box.min_x; x <= tri_bounding_box.max_x; ++x) {
                    std::cout << framebuffer[x - tri_bounding_box.min_x][y - tri_bounding_box.min_y]; // print the framebuffer if a point exists, add offset
                }
                std::cout << "" << std::endl; // move to the next row
            }
        }
};

class Operations3D {
    public:
        static Point project3DPoint(Point3D p, double focus) {
            Point projected = Point(p.x * focus / p.z, p.y * focus / p.z);
            return projected;
        }
        static std::vector<Triangle> project3DTriangleArr(std::vector<Triangle3D> Triangles3D, int length, double focus) {
            std::vector<Triangle> triangles;
            for (size_t k = 0; k < length; k++) {
                    // Trim triangles that go behind
                    if (Triangles3D[k].A.z < 0 || Triangles3D[k].B.z < 0 || Triangles3D[k].C.z < 0) {
                        continue;
                    }
                    Point pointA = project3DPoint(Triangles3D[k].A, focus);
                    Point pointB = project3DPoint(Triangles3D[k].B, focus);
                    Point pointC = project3DPoint(Triangles3D[k].C, focus);
                    Triangle object_to_push_back = Triangle(pointA, pointB, pointC);
                    object_to_push_back.material = Triangles3D[k].material;
                    triangles.push_back(object_to_push_back);
                }
            return triangles;
        }

        static std::vector<Triangle3D> transform3DTriangleArr(std::vector<Triangle3D> Triangles3D, double length) {
            double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count();
            seconds *= 1.5;
            double s = sin(seconds);
            double c = cos(seconds);
            std::vector<Triangle3D> output_triangles;
            Matrix rotatez = Matrix(c, 0, s, 0, 1, 0, -s, 0, c);
            Matrix rotatex = Matrix(1, 0, 0, 0, c, -s, 0, s, c);
            for (size_t k = 0; k < length; k++) {
                Point3D point1 = rotatez.multiply_with_vector(Triangles3D[k].A);
                Point3D point2 = rotatez.multiply_with_vector(Triangles3D[k].B);
                Point3D point3 = rotatez.multiply_with_vector(Triangles3D[k].C);
                point1 = rotatex.multiply_with_vector(point1);
                point2 = rotatex.multiply_with_vector(point2);
                point3 = rotatex.multiply_with_vector(point3);
                // translate them back so they aren't right in front of the camera
                point1.z += transl_distance;
                point2.z += transl_distance;
                point3.z += transl_distance;
                output_triangles.push_back(Triangle3D(point1, point2, point3));
            }
            return output_triangles;
        }
        static std::vector<Triangle3D> sortAndMaterializeTriangleArr(std::vector<Triangle3D> Triangles3D, double length) {
            std::string shades = "-~+*#@$";
            for (size_t k = 0; k < length; k++) {
                // compute B - A and C - A
                Point3D vec1 = Triangles3D[k].B.add(Triangles3D[k].A.negate()).normalize();
                Point3D vec2 = Triangles3D[k].C.add(Triangles3D[k].A.negate()).normalize();
                // compute normalized cross product aka the direction the triangle is facing
                Point3D normal = vec1.cross(vec2).normalize();

                // now this is the hardcoded (for now) light direction vector (1, 1, 1) normalized
                Point3D direction = Point3D(0, 0, -1);
                // compute dot
                double dotprod = direction.dot(normal);
                // scale it up
                dotprod = std::max(dotprod, 0.0);
                dotprod *= 7;
                // now pick the shade value, values close to 1 will have lighter shade
                Triangles3D[k].material = shades[(int) dotprod];
            }

            std::sort(Triangles3D.begin(), Triangles3D.end(), compareAvgZVal);
            return Triangles3D;
        }
        static bool compareAvgZVal(Triangle3D w, Triangle3D u) {
            double avgwZ = (w.A.z + w.B.z + w.C.z) / 3;
            double avguZ = (u.A.z + u.B.z + u.C.z) / 3;
            return avgwZ > avguZ;
        }
};

class MeshLoader {
    public:
    static std::vector<Triangle3D> loadFromObjFile(std::string file_path) {
        std::vector <std::vector<double>> output_vertices; // list of list of numbers (list of vertices)
        std::vector <std::vector<int>> output_faces;
        std::vector<double> temp(3);
        std::vector<int> temp_faces(3);

        std::ifstream file = std::ifstream(file_path);
        if (!file) {
            std::cerr << "Failed to open file.\n";
            exit(1);
        }

        std::string line;
        while (std::getline(file, line)) {  // read line by line
            int commentedOut = line.find("#");
            if (commentedOut != std::string::npos) {
                line = line.erase(commentedOut);
            }
            // if the line is skipped
            if (line.empty()) {
                continue;
            }
            if (line.at(0) == 'v' && line.at(1) == ' ') {
                line.replace(0, 2, "");
                std::vector<std::string> points = split(line, ' ');
                // remove any extraneous entries
                points.erase(std::remove(points.begin(), points.end(), ""), points.end());

                temp =  {std::stod(points[0]), std::stod(points[1]), std::stod(points[2])};
                output_vertices.push_back(temp);
            }
            if (line.at(0) == 'f') {
                line.replace(0, 2, "");
                std::vector<std::string> points = split(line, ' ');
                // remove any extranenous entries
                points.erase(std::remove(points.begin(), points.end(), ""), points.end());
                
                // ignore anything after the slash
                for (size_t k = 0; k < points.size(); k++) {
                    points[k] = split(points[k], '/')[0];
                }

                temp_faces =  {std::stoi(points[0]), std::stoi(points[1]), std::stoi(points[2])};
                output_faces.push_back(temp_faces);
            }
        }
        std::vector<Triangle3D> output_triangles;
        for (size_t k = 0; k < output_faces.size(); k++) {
            // convert it to an array of triangles object
            // get the face, and look up the vertex for each one
            // obj files are indexed from 1 so subtract one
            std::vector<double> p1vec = output_vertices[output_faces[k][0] - 1];
            Point3D p1 = Point3D(p1vec[0], p1vec[1], p1vec[2]);

            std::vector<double> p2vec = output_vertices[output_faces[k][1] - 1];
            Point3D p2 = Point3D(p2vec[0], p2vec[1], p2vec[2]);

            std::vector<double> p3vec = output_vertices[output_faces[k][2] - 1];
            Point3D p3 = Point3D(p3vec[0], p3vec[1], p3vec[2]);
            output_triangles.push_back(Triangle3D(p1, p2, p3));
        }
        return output_triangles;
    }
    // helper function to split a string
    static std::vector<std::string> split(const std::string& s, char delimiter) {
        std::vector<std::string> tokens;
        size_t start = 0;
        size_t end = s.find(delimiter);

        while (end != std::string::npos) {
            tokens.push_back(s.substr(start, end - start));
            start = end + 1;
            end = s.find(delimiter, start);
        }
        tokens.push_back(s.substr(start)); // last token
        return tokens;
    }

};
int main() {
    std::string path;
    std::string zoomstr;
    std::string widthminstr;
    std::string heightminstr;
    std::string widthmaxstr;
    std::string heightmaxstr;
    std::string transl_distance_str;
    double zoom;
    int width_min;
    int height_min;
    int width_max;
    int height_max;
    std::cout << "What would you like to render (path to .obj file): " << std::endl;
    std::cin >> path;
    std::cout << "Zoom value?" << std::endl;
    std::cin >> zoomstr;
    std::cout << "Enter screen width min: " << std::endl;
    std::cin >> widthminstr;
    std::cout << "Enter screen height min: " << std::endl;
    std::cin >> heightminstr;
    std::cout << "Enter screen width max: " << std::endl;
    std::cin >> widthmaxstr;
    std::cout << "Enter screen height max: " << std::endl;
    std::cin >> heightmaxstr;
    std::cout << "Enter translation distance (5 works for most small 3d objects): " << std::endl;
    std::cin >> transl_distance_str;

    height_max = std::stoi(heightmaxstr);
    width_max = std::stoi(widthmaxstr);
    height_min = std::stoi(heightminstr);
    width_min = std::stoi(widthminstr);
    zoom = std::stod(zoomstr);
    transl_distance = std::stod(transl_distance_str);

    std::vector<Triangle3D> triangle_to_render = MeshLoader::loadFromObjFile(path);
    std::srand(std::time(0));

    std::vector<Triangle3D> new_triangle(triangle_to_render.size());
    Board b1 = Board(width_min, height_min, width_max, height_max);
    std::vector<Triangle3D> transformed_triangles;
    std::vector<Triangle3D> materialized_triangles;
    std::vector<Triangle> projected_triangles;
    while (1) {
        std::this_thread::sleep_for(std::chrono::milliseconds(16));
        for (size_t k = 0; k < triangle_to_render.size(); k++) {
            new_triangle[k] = triangle_to_render[k];
        }
        transformed_triangles = Operations3D::transform3DTriangleArr(new_triangle, new_triangle.size());
        materialized_triangles = Operations3D::sortAndMaterializeTriangleArr(transformed_triangles, transformed_triangles.size());
        projected_triangles = Operations3D::project3DTriangleArr(materialized_triangles, materialized_triangles.size(), zoom);
        b1.print_board(projected_triangles, projected_triangles.size());
    }
    return 0;
}
