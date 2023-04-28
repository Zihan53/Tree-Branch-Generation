#include <iostream>
#include <vector>
#include <ctime>
#include <Eigen/Core>
#include <random>
#include <igl/segment_segment_intersect.h>

using namespace Eigen;
//0.590162
double RESIZE_RATIO = 0.9;
int RESIZE_MAX_TIME = 8;
int GEN_ROOT_NUM = 3;
int SMALLEST_BRANCH_SIZE = 30;
int MONTE_CARLO_TIME = 30000;
std::vector<double> left_angle = {M_PI/4.0, M_PI/12.0};
std::vector<double> middle_angle = {-M_PI/18.0, M_PI/18.0};
std::vector<double> right_angle = {-M_PI/12.0, -M_PI/4.0};
std::vector<double> ratio_list = {0.8, 0.9, 0.8};
class Root {
public:
    Vector3d point;
    Vector3d direction;
    double height;
    int lmf; //0:left, 1:middle, 2:right
    Root(const Vector3d& p, const Vector3d& d, double h, int i): point(p), direction(d), height(h), lmf(i) {};
};

class Branch {
public:
    double height;
    double width;
    Root root;

    Vector3d left_top;
    Vector3d right_top;
    Vector3d left_bottom;
    Vector3d right_bottom;

    Branch(double h, double w, const Root& r): height(h), width(w), root(r) {
        Vector3d base_vector(-root.direction.y()*width, root.direction.x()*width, 0);
        Vector3d height_vector = root.direction*height;
        left_bottom = root.point+base_vector/2;
        right_bottom = root.point-base_vector/2;
        left_top = root.point+base_vector/2+height_vector;
        right_top = root.point-base_vector/2+height_vector;
    }

    std::vector<Root> generate_roots() {
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<double> leftdis(left_angle[0], left_angle[1]);
        std::uniform_real_distribution<double> middledis(middle_angle[0], middle_angle[1]);
        std::uniform_real_distribution<double> rightdis(right_angle[0], right_angle[1]);
        std::vector<Root> roots;
        Vector3d direction = root.direction;
        Vector3d point = root.point;
        Vector3d left_direction(std::cos(leftdis(gen))*direction.x()-std::sin(leftdis(gen))*direction.y(), std::sin(leftdis(gen))*direction.x()+std::cos(leftdis(gen))*direction.y(), 0);
        Vector3d right_direction(std::cos(rightdis(gen))*direction.x()-std::sin(rightdis(gen))*direction.y(), std::sin(rightdis(gen))*direction.x()+std::cos(rightdis(gen))*direction.y(), 0);
        Vector3d middle_direction(std::cos(middledis(gen))*direction.x()-std::sin(middledis(gen))*direction.y(), std::sin(middledis(gen))*direction.x()+std::cos(middledis(gen))*direction.y(), 0);
        Root left_root(point+height*direction, left_direction, height*ratio_list[0], 0);
        Root right_root(point+height*direction, right_direction, height*ratio_list[2], 2);
        Root middle_root(point+height*direction, middle_direction, height*ratio_list[1], 1);
        if (left_root.height > SMALLEST_BRANCH_SIZE*ratio_list[left_root.lmf]) {
            roots.push_back(left_root);
        }
        if (left_root.height > SMALLEST_BRANCH_SIZE*ratio_list[right_root.lmf]) {
            roots.push_back(right_root);
        }
        if (left_root.height > SMALLEST_BRANCH_SIZE*ratio_list[middle_root.lmf]) {
            roots.push_back(middle_root);
        }
        return roots;
    }

    bool check_collision_branch(const Branch& b) {
        double t = 0, u = 0;
        bool intersects = false;
        intersects = intersects || igl::segment_segment_intersect(left_bottom, left_top - left_bottom, b.left_bottom, b.left_top - b.left_bottom, t, u);
        intersects = intersects || igl::segment_segment_intersect(left_bottom, left_top - left_bottom, b.left_top, b.right_top - b.left_top, t, u);
        intersects = intersects || igl::segment_segment_intersect(left_bottom, left_top - left_bottom, b.right_top, b.right_bottom - b.right_top, t, u);
        //intersects = intersects || igl::segment_segment_intersect(left_bottom, left_top - left_bottom, b.right_bottom, b.left_bottom - b.right_bottom, t, u);

        intersects = intersects || igl::segment_segment_intersect(left_top, right_top - left_top, b.left_bottom, b.left_top - b.left_bottom, t, u);
        intersects = intersects || igl::segment_segment_intersect(left_top, right_top - left_top, b.left_top, b.right_top - b.left_top, t, u);
        intersects = intersects || igl::segment_segment_intersect(left_top, right_top - left_top, b.right_top, b.right_bottom - b.right_top, t, u);
        //intersects = intersects || igl::segment_segment_intersect(left_top, right_top - left_top, b.right_bottom, b.left_bottom - b.right_bottom, t, u);

        intersects = intersects || igl::segment_segment_intersect(right_top, right_bottom - right_top, b.left_bottom, b.left_top - b.left_bottom, t, u);
        intersects = intersects || igl::segment_segment_intersect(right_top, right_bottom - right_top, b.left_top, b.right_top - b.left_top, t, u);
        intersects = intersects || igl::segment_segment_intersect(right_top, right_bottom - right_top, b.right_top, b.right_bottom - b.right_top, t, u);
        //intersects = intersects || igl::segment_segment_intersect(right_top, right_bottom - right_top, b.right_bottom, b.left_bottom - b.right_bottom, t, u);

        // intersects = intersects || igl::segment_segment_intersect(right_bottom, left_bottom - right_bottom, b.left_bottom, b.left_top - b.left_bottom, t, u);
        // intersects = intersects || igl::segment_segment_intersect(right_bottom, left_bottom - right_bottom, b.left_top, b.right_top - b.left_top, t, u);
        // intersects = intersects || igl::segment_segment_intersect(right_bottom, left_bottom - right_bottom, b.right_top, b.right_bottom - b.right_top, t, u);
        // intersects = intersects || igl::segment_segment_intersect(right_bottom, left_bottom - right_bottom, b.right_bottom, b.left_bottom - b.right_bottom, t, u);
        return intersects;
    }

    bool check_collision_canopy(const Vector3d& start, const Vector3d& end) {
        double t = 0, u = 0;
        bool intersects = false;
        intersects = intersects || igl::segment_segment_intersect(left_bottom, left_top - left_bottom, start, end - start, t, u);
        intersects = intersects || igl::segment_segment_intersect(left_top, right_top - left_top, start, end - start, t, u);
        intersects = intersects || igl::segment_segment_intersect(right_top, right_bottom - right_top, start, end - start, t, u);
        //intersects = intersects || igl::segment_segment_intersect(right_bottom, left_bottom - right_bottom, start, end - start, t, u);
        return intersects;
    }

    void resize(double ratio) {
        height *= ratio;
        width *= ratio;
        Vector3d height_vector = root.direction*height;
        Vector3d base_vector = (left_bottom-right_bottom)*ratio;
        left_bottom = root.point+base_vector/2;
        right_bottom = root.point-base_vector/2;
        left_top = root.point+base_vector/2+height_vector;
        right_top = root.point-base_vector/2+height_vector;
    }
};

class Tree {
public:
    std::vector<Branch> branches;
    std::vector<Vector3d> canopy;
    Tree(double mh, double mw, std::vector<Vector3d>& c): canopy(c) {
        Vector3d direction(0, 1, 0);
        Vector3d start_point(-100, -150, 0);
        Root start_root(start_point, direction, mh, 1);
        Branch branch(mh, mw, start_root);
        branches.push_back(branch);

        std::vector<Root> generated_roots = branch.generate_roots();
        roots.insert(roots.end(), generated_roots.begin(), generated_roots.end());
    }

    void generate() {
        while (roots.size() > 0) {
            std::cout << roots.size() << std::endl;
            // Randomly selected a root and check branch size is large enough
            int index = rand() % roots.size();
            Root curr_root = roots[index];
            roots.erase(roots.begin() + index);
            if (curr_root.height < SMALLEST_BRANCH_SIZE*ratio_list[curr_root.lmf])
                continue;

            // Find the branch with proper size and check branch size is large enough
            Branch new_branch = find_branch(curr_root.height, curr_root);
            if (new_branch.height >= SMALLEST_BRANCH_SIZE*ratio_list[curr_root.lmf]) {
                branches.push_back(new_branch);
                std::vector<Root> curr_roots = new_branch.generate_roots();
                roots.insert(roots.end(), curr_roots.begin(), curr_roots.end());
            }
        }
    }

    void generate_one_branch() {
        int curr_num = branches.size();
        std::cout << roots.size() << std::endl;
        while (branches.size() == curr_num && roots.size() > 0) {
            int index = rand() % roots.size();
            Root curr_root = roots[index];
            roots.erase(roots.begin() + index);
            if (curr_root.height < SMALLEST_BRANCH_SIZE*ratio_list[curr_root.lmf])
                continue;

            // Find the branch with proper size and check branch size is large enough
            Branch new_branch = find_branch(curr_root.height, curr_root);
            if (new_branch.height >= SMALLEST_BRANCH_SIZE*ratio_list[curr_root.lmf]) {
                branches.push_back(new_branch);
            }
        }
    }

    void simulated_annealing() {
        double curr_rate = rate();

        std::vector<Branch> best_tree = branches;
        double best_rate = rate();

        int i = 0;
        while (roots.size() > 0) {
            std::cout<<"\n";
            generate_one_branch(); //Give a purturbing of current tree
            double new_rate = rate();
            double delta = new_rate - curr_rate;
            std::cout << new_rate << std::endl;
            std::cout << best_rate << std::endl;

            std::random_device rd;
            std::mt19937 gen(rd());
            std::uniform_real_distribution<double> uniform_dist(0.0, 1.0);
            if (delta < 0 || uniform_dist(gen) < 0.5) {
                curr_rate = new_rate;
                if (new_rate < best_rate) {
                    best_tree = branches;
                    best_rate = new_rate;
                }
                std::vector<Root> curr_roots = branches.back().generate_roots();
                roots.insert(roots.end(), curr_roots.begin(), curr_roots.end());
            } else {
                //roots.push_back(branches.back().root);
                branches.pop_back();
            }
            i++;
        }
        branches = best_tree;
        std::cout << branches.size() << std::endl;
    } 

    double rate() {
        double a = area();
        double ar = avg_ratio();
        std::cout << a << " " << ar << std::endl;
        return 0.7 * a + ar * 0.3;
    }

private:
    std::vector<Root> roots; //Generated roots

    double area() {
        // Find the bounding square
        auto xminmax = std::minmax_element(canopy.begin(), canopy.end(), [](const Vector3d& a, const Vector3d& b) { return a.x() < b.x(); });
        auto yminmax = std::minmax_element(canopy.begin(), canopy.end(), [](const Vector3d& a, const Vector3d& b) { return a.y() < b.y(); });
        double xmin, xmax, ymin, ymax;
        xmin = (*xminmax.first)[0];
        xmax = (*xminmax.second)[0];
        ymin = (*yminmax.first)[1];
        ymax = (*yminmax.second)[1];

        // Monte Carlo method
        int canopy_count = 0;
        int branch_count = 0;
        for (auto i=0; i<MONTE_CARLO_TIME; i++) {
            std::random_device rd;
            std::mt19937 gen(rd());
            std::uniform_real_distribution<double> xdis(xmin, xmax);
            std::uniform_real_distribution<double> ydis(ymin, ymax);
            Vector3d point(xdis(gen), ydis(gen), 0);
            if(point_inside_shape(canopy, point)) {
                canopy_count++;
                for (auto i=1; i<branches.size(); i++) {
                    std::vector<Vector3d> branchList;
                    branchList.push_back(branches[i].left_bottom);
                    branchList.push_back(branches[i].right_bottom);
                    branchList.push_back(branches[i].right_top);
                    branchList.push_back(branches[i].left_top);
                    if(point_inside_shape(branchList, point)) {
                        branch_count++;
                        break;
                    }
                }
            }
        }
        return 1.0 * (canopy_count-branch_count) / canopy_count;
    }

    double avg_ratio() {
        double ratio = 1;
        if (branches.size() == 1)
            return ratio;
        for (auto i=1; i<branches.size(); i++) {
            double curr_ratio = 1 - branches[i].height / (branches[i].root.height);
            // if (curr_ratio>0.5)
            //     curr_ratio *= 5;
            ratio += curr_ratio;
            //std::cout<<"hhh"<<branches[i].height<<" "<<branches[i].root.height<<std::endl;
        }
        return ratio / (branches.size() - 1);
    }

    Branch find_branch(double max_height, Root& curr_root) {
        bool intersects;
        Branch new_branch(max_height, max_height/7, curr_root);
        for (auto j=0; j<RESIZE_MAX_TIME; j++) {
            intersects = false;
            // Check intersection with branches
            for (const Branch& b: branches) {
                if (b.root.point != new_branch.root.point && b.root.point+b.height*b.root.direction != new_branch.root.point) {
                    intersects = new_branch.check_collision_branch(b);
                }
                if (intersects)
                    break;
            }
            // Check intersection with canopy
            for (auto i=0; i<canopy.size(); i++) {
                if (intersects)
                    break;
                intersects = new_branch.check_collision_canopy(canopy[i], canopy[(i+1) % canopy.size()]);
            }
            // Resize
            if (intersects) {
               new_branch.resize(RESIZE_RATIO);
            } else {
               break;
            }
            if (new_branch.height<SMALLEST_BRANCH_SIZE*ratio_list[curr_root.lmf]) {
                intersects = true;
                break;
            }
        }
        if (intersects) {
            new_branch.resize(0);
        }
        return new_branch;
    }

    bool point_inside_shape(std::vector<Vector3d>& shape, Vector3d point) {
        int num_intersections = 0;
        for (int i = 0; i < shape.size(); i++) {
            Vector3d p1 = shape[i];
            Vector3d p2 = shape[(i+1) % shape.size()];
            if (point_on_edge(point, p1, p2)) {
                return true;
            }
            if (edge_intersects_ray(point, p1, p2)) {
                num_intersections++;
            }
        }
        return num_intersections % 2 == 1;
    }

    bool point_on_edge(Vector3d point, Vector3d edge_start, Vector3d edge_end) {
        if (point == edge_start || point == edge_end) {
            return true;
        }
        double dx = edge_end[0] - edge_start[0];
        double dy = edge_end[1] - edge_start[1];
        if (dx == 0) {
            return point[0] == edge_start[0] && edge_start[1] <= point[1] && point[1] <= edge_end[1];
        } else if (dy == 0) {
            return point[1] == edge_start[1] && edge_start[0] <= point[0] && point[0] <= edge_end[0];
        } else {
            double t1 = (point[0] - edge_start[0]) / dx;
            double t2 = (point[1] - edge_start[1]) / dy;
            return abs(t1 - t2) < 0.00001 && 0 <= t1 && t1 <= 1 && 0 <= t2 && t2 <= 1;
        }
    }

    bool edge_intersects_ray(Vector3d point, Vector3d edge_start, Vector3d edge_end) {
        if (point[1] < std::min(edge_start[1], edge_end[1]) || point[1] > std::max(edge_start[1], edge_end[1]) || point[1] == edge_end[1]) {
            return false;
        } 
        if (point[0] > std::max(edge_start[0], edge_end[0])) {
            return false;
        }
        if (point[0] < std::min(edge_start[0], edge_end[0])) {
            return true;
        }
        double slope = (edge_end[0] - edge_start[0]) / (edge_end[1] - edge_start[1]);
        double x_intercept = (point[1] - edge_start[1]) * slope + edge_start[0];
        return x_intercept >= point[1];
    }
};