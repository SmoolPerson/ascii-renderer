#include <iostream>
#include <ctime>
#include <cmath>
#include <thread>
#include <vector>
#include <algorithm>
#include <fstream>
#include <SFML/Graphics.hpp>
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
    bool operator==(const Point& other) const {
        return x == other.x && y == other.y;
    }
};

class Triangle {
    public:
        Point A;
        Point B;
        Point C;
        int shade = 0;
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
            if (fabs(s) < 0.01 || fabs(t) < 0.01 || fabs(s + t - 1) < 0.01) {
                return true;
            }
            return false;
        }
        bool operator==(const Triangle& other) const {
            return (A == other.A && B == other.B && C == other.C) ||
           (A == other.B && B == other.C && C == other.A) ||
           (A == other.C && B == other.A && C == other.B);
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
        int shade = 0;

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
class Screen {
    public:
        BoundingBox tri_bounding_box;
        sf::RenderWindow window;
        sf::Image image;
        sf::Texture texture;
        uint width;
        uint height;
        Screen(double min_x, double min_y, double max_x, double max_y) {
            tri_bounding_box.min_x = min_x;
            tri_bounding_box.max_x = max_x;
            tri_bounding_box.max_y = max_y;
            tri_bounding_box.min_y = min_y;
            width = abs(max_x - min_x);
            height = abs(max_y - min_y);
            window = sf::RenderWindow(sf::VideoMode({width, height}), "SFML window");
            window.setFramerateLimit(60);
        }

        void print_board(std::vector<Triangle> triangle, int length) {
            image = sf::Image({width, height}, sf::Color::Black);
            int total_iters = 0;
            auto startTime = std::chrono::steady_clock::now();
            for (size_t l = 0; l < triangle.size(); l++) {
                int temp_iters = 0;
                BoundingBox temp; // temporary bounding box for the bounding box of the current triangle so we loop over less points

                // compute bounding box then loop over that
                temp.max_x = std::max(triangle[l].A.x, std::max(triangle[l].B.x, triangle[l].C.x));
                temp.max_y = std::max(triangle[l].A.y, std::max(triangle[l].B.y, triangle[l].C.y));
                temp.min_x = std::min(triangle[l].A.x, std::min(triangle[l].B.x, triangle[l].C.x));
                temp.min_y = std::min(triangle[l].A.y, std::min(triangle[l].B.y, triangle[l].C.y));

                // clamp it so it doesnt go outside the screen width/height
                temp.min_x = std::max(temp.min_x, tri_bounding_box.min_x);
                temp.max_x = std::min(temp.max_x, tri_bounding_box.max_x);
                temp.min_y = std::max(temp.min_y, tri_bounding_box.min_y);
                temp.max_y = std::min(temp.max_y, tri_bounding_box.max_y);
                for (int x = temp.min_x; x <= temp.max_x; x++) {
                    for (int y = temp.min_y; y <= temp.max_y; y++) {
                        total_iters++;
                        temp_iters++;
                        if (triangle[l].point_in_triangle(Point(x, y))) {
                            sf::Color col = image.getPixel({(uint)(x - tri_bounding_box.min_x), (uint) (y - tri_bounding_box.min_y)});
                            if (col == sf::Color::Black) {
                                image.setPixel({(uint)(x - tri_bounding_box.min_x), (uint) (y - tri_bounding_box.min_y)}, sf::Color(triangle[l].shade, triangle[l].shade, triangle[l].shade));
                            }
                        }
                    }
                }
            }
            auto endTime = std::chrono::steady_clock::now();
            auto duration = std::chrono::duration<double>(endTime - startTime);
            while (const std::optional event = window.pollEvent())
            {
            // Close window: exit
            if (event->is<sf::Event::Closed>())
                window.close();
            }

            window.clear();
            texture.loadFromImage(image);
            sf::Sprite sprite(texture);
            window.draw(sprite);
            window.display();
        }
};


class Player {
    public:
        Point3D current_pos;
        double z_angle = 0.1;
        double x_angle = 0.1;
        Player(Point3D pos) {
            current_pos = pos;
        }
        Player() {
            current_pos = Point3D(0, 0, 5);
        }
        Point3D rotateAmtZ(Point3D rotationAmt) {
            double s = sin(-z_angle);
            double c = cos(-z_angle);
            Matrix rotatez = Matrix(c, 0, s, 0, 1, 0, -s, 0, c); // z rotation matrix
            Point3D rotated = rotatez.multiply_with_vector(rotationAmt);
            return rotated;
        }
        void handle_movement() {
            if (sf::Keyboard::isKeyPressed(sf::Keyboard::Key::W)) {
                current_pos = current_pos.add(rotateAmtZ(Point3D(0, 0, -0.06))); // add the move amount rotated by the angle
            }
            if (sf::Keyboard::isKeyPressed(sf::Keyboard::Key::S)) {
                current_pos = current_pos.add(rotateAmtZ(Point3D(0, 0, 0.06)));
            }
            if (sf::Keyboard::isKeyPressed(sf::Keyboard::Key::A)) {
                current_pos = current_pos.add(rotateAmtZ(Point3D(0.06, 0, 0)));
            }
            if (sf::Keyboard::isKeyPressed(sf::Keyboard::Key::D)) {
                current_pos = current_pos.add(rotateAmtZ(Point3D(-0.06, 0, 0)));
            }
            if (sf::Keyboard::isKeyPressed(sf::Keyboard::Key::Space)) {
                current_pos.y -= 0.06;
            }
            if (sf::Keyboard::isKeyPressed(sf::Keyboard::Key::LShift)) {
                current_pos.y += 0.06;
            }
            if (sf::Keyboard::isKeyPressed(sf::Keyboard::Key::Left)) {
                z_angle += 0.04;
            }
            if (sf::Keyboard::isKeyPressed(sf::Keyboard::Key::Right)) {
                z_angle -= 0.04;
            }
            if (sf::Keyboard::isKeyPressed(sf::Keyboard::Key::Up)) {
                x_angle += 0.04;
            }
            if (sf::Keyboard::isKeyPressed(sf::Keyboard::Key::Down)) {
                x_angle -= 0.04;
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
                    // flip y to make the image be right-side up
                    Triangles3D[k].A.y = -Triangles3D[k].A.y;
                    Triangles3D[k].B.y = -Triangles3D[k].B.y;
                    Triangles3D[k].C.y = -Triangles3D[k].C.y;
                    Point pointA = project3DPoint(Triangles3D[k].A, focus);
                    Point pointB = project3DPoint(Triangles3D[k].B, focus);
                    Point pointC = project3DPoint(Triangles3D[k].C, focus);
                    Triangle object_to_push_back = Triangle(pointA, pointB, pointC);
                    object_to_push_back.shade = Triangles3D[k].shade;
                    triangles.push_back(object_to_push_back);
                }
            return triangles;
        }

        static std::vector<Triangle3D> transform3DTriangleArr(std::vector<Triangle3D> Triangles3D, double length, Player player) {
            double z_angle = player.z_angle;
            double x_angle = player.x_angle;
            double s = sin(z_angle);
            double c = cos(z_angle);
            double sx = sin(x_angle);
            double cx = cos(x_angle);
            std::vector<Triangle3D> output_triangles;
            Matrix rotatez = Matrix(c, 0, s, 0, 1, 0, -s, 0, c);
            Matrix rotatex = Matrix(1, 0, 0, 0, cx, -sx, 0, sx, cx);
            for (size_t k = 0; k < length; k++) {
                Point3D point1 = Triangles3D[k].A;
                Point3D point2 = Triangles3D[k].B;
                Point3D point3 = Triangles3D[k].C;
                // translate them based on the players position
                point1.z += player.current_pos.z;
                point2.z += player.current_pos.z;
                point3.z += player.current_pos.z;
                point1.x += player.current_pos.x;
                point2.x += player.current_pos.x;
                point3.x += player.current_pos.x;
                point1.y += player.current_pos.y;
                point2.y += player.current_pos.y;
                point3.y += player.current_pos.y;
                point1 = rotatez.multiply_with_vector(point1);
                point2 = rotatez.multiply_with_vector(point2);
                point3 = rotatez.multiply_with_vector(point3);
                point1 = rotatex.multiply_with_vector(point1);
                point2 = rotatex.multiply_with_vector(point2);
                point3 = rotatex.multiply_with_vector(point3);
                output_triangles.push_back(Triangle3D(point1, point2, point3));
            }
            return output_triangles;
        }
        static Point3D getNormal(Triangle3D tri) {
            // compute B - A and C - A to get two vectors that align the triangle (normalize since we only care about direction)
            Point3D vec1 = tri.B.add(tri.A.negate()).normalize();
            Point3D vec2 = tri.C.add(tri.A.negate()).normalize();
            Point3D normal = vec1.cross(vec2).normalize(); // do cross product to get the normal vector
            return normal;
        }
        static std::vector<Triangle3D> sortAndMaterializeTriangleArr(std::vector<Triangle3D> Triangles3D, double length, Player player) {
            std::vector<Triangle3D> newTriangles3D;
            for (size_t k = 0; k < length; k++) {
                Point3D normal = getNormal(Triangles3D[k]);
                // now this is the hardcoded (for now) light direction vector (0, 0, -1) normalized
                Point3D direction = Point3D(0, 0, -1);
                // compute dot
                double dotprod = direction.dot(normal);
                // rotate current position so it matches the angle
                double c = cos(player.z_angle);
                double s = sin(player.z_angle);
                Matrix rotatez = Matrix(c, 0, s, 0, 1, 0, -s, 0, c);
                // scale it up, max of 0.3 so there is ambient lighting
                dotprod = std::max(dotprod, 0.3);
                // now pick the shade value, values close to 1 will have lighter shade
                newTriangles3D.push_back(Triangles3D[k]);
                newTriangles3D[newTriangles3D.size() - 1].shade = (int) (dotprod*256); // accessing latest element, and shade is from 0-255 (rgb)
            }
            std::sort(newTriangles3D.begin(), newTriangles3D.end(), compareAvgZVal);
            return newTriangles3D;
        }
        static bool compareAvgZVal(Triangle3D w, Triangle3D u) {
            double avgwZ = (w.A.z + w.B.z + w.C.z);
            double avguZ = (u.A.z + u.B.z + u.C.z);
            return avgwZ < avguZ;
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
    std::string widthminstr;
    std::string heightminstr;
    std::string widthmaxstr;
    std::string heightmaxstr;
    int width_min;
    int height_min;
    int width_max;
    int height_max;
    std::cout << "What would you like to render (path to .obj file): " << std::endl;
    std::cin >> path;
    std::cout << "Enter screen width min: " << std::endl;
    std::cin >> widthminstr;
    std::cout << "Enter screen height min: " << std::endl;
    std::cin >> heightminstr;
    std::cout << "Enter screen width max: " << std::endl;
    std::cin >> widthmaxstr;
    std::cout << "Enter screen height max: " << std::endl;
    std::cin >> heightmaxstr;

    height_max = std::stoi(heightmaxstr);
    width_max = std::stoi(widthmaxstr);
    height_min = std::stoi(heightminstr);
    width_min = std::stoi(widthminstr);

    std::vector<Triangle3D> triangle_to_render = MeshLoader::loadFromObjFile(path);
    Player player1 = Player();
    Screen b1 = Screen(width_min, height_min, width_max, height_max);
    std::vector<Triangle3D> transformed_triangles;
    std::vector<Triangle3D> materialized_triangles;
    std::vector<Triangle> projected_triangles;
    while (1) {
        player1.handle_movement();
        transformed_triangles = Operations3D::transform3DTriangleArr(triangle_to_render, triangle_to_render.size(), player1);
        materialized_triangles = Operations3D::sortAndMaterializeTriangleArr(transformed_triangles, transformed_triangles.size(), player1);
        projected_triangles = Operations3D::project3DTriangleArr(materialized_triangles, materialized_triangles.size(), 500);\
        b1.print_board(projected_triangles, projected_triangles.size());
        if (!b1.window.isOpen()) {
            return 1;
        }
    }
    return 0;
}
