/*
    A multi-threaded ray tracing demo using GLFW for windowing.

    Brandon Luk
*/

#define _USE_MATH_DEFINES

#include "color.hpp"
#include "ray_math.hpp"
#include "entity.hpp"
#include "vec3.hpp"

#include <GLFW/glfw3.h>

#include <algorithm>
#include <thread>
#include <chrono>
#include <cmath>
#include <iostream>
#include <limits>
#include <string>
#include <memory>
#include <sstream>
#include <utility>
#include <vector>


std::vector<std::thread> thread_pool;

std::string WINDOW_TITLE = "Raytracer";

constexpr int VIEW_FRAME_WIDTH = 1280;
constexpr int VIEW_FRAME_HEIGHT = 960;

constexpr int COLOR_DEPTH = 3; // 3 since we are using BGR color representation
constexpr int COLOR_VALUE_MAX = 255;

constexpr int REFLECTIVE_DEPTH_MAX = 5; // Maximum number of times we will bouce a ray off reflective surfaces

constexpr Color BACKGROUND_COLOR{200, 155, 55};

constexpr int PIXEL_DATA_SIZE = VIEW_FRAME_WIDTH * VIEW_FRAME_HEIGHT * COLOR_DEPTH;
GLubyte* pixel_data = new GLubyte[PIXEL_DATA_SIZE];

const double FOV = 30.0 * M_PI / 180.0;

struct Light{
    Vec3<double> origin;
    double intensity;
};

struct{
    double offset = 5.0;
    Vec3<double> normal = {0.0, 0.0 , 1.0};
} view_frame;


Vec3<double> camera = {0.0, 0.0, 0.0};
constexpr double CAMERA_MOVE_SPEED = 10.0;
constexpr double CAMERA_ROTATE_SPEED = 0.5;


Vec3<double> roll_component = {0.0, 1.0, 0.0}; // Defines which direction we will be considering as "up"

std::vector<std::unique_ptr<Entity> > entities;

std::vector<Light> lights{   {{0.0, 10.0, 10.0}, 1.3}};


// Keyboard controls
/**************************************************************************************************************/

// Key press flags
bool W_KEY_PRESSED = false;
bool A_KEY_PRESSED = false;
bool S_KEY_PRESSED = false;
bool D_KEY_PRESSED = false;

bool UP_KEY_PRESSED = false;
bool LEFT_KEY_PRESSED = false;
bool RIGHT_KEY_PRESSED = false;
bool DOWN_KEY_PRESSED = false;

void MoveCameraForward(double time_delta)
{
    camera = camera + (time_delta * CAMERA_MOVE_SPEED * view_frame.normal);
}

void MoveCameraLeft(double time_delta)
{
    Vec3<double> rotated_view_frame = view_frame.normal;
    rotated_view_frame = rotated_view_frame.Rotate_y(M_PI_2);

    camera = camera - (time_delta * CAMERA_MOVE_SPEED * rotated_view_frame);
}

void MoveCameraRight(double time_delta)
{
    Vec3<double> rotated_view_frame = view_frame.normal;
    rotated_view_frame = rotated_view_frame.Rotate_y(M_PI_2);

    camera = camera + (time_delta * CAMERA_MOVE_SPEED * rotated_view_frame);
}

void MoveCameraBackward(double time_delta)
{
    camera = camera - (time_delta * CAMERA_MOVE_SPEED * view_frame.normal);
}

void RotateCameraUpwards(double time_delta)
{
    view_frame.normal = view_frame.normal.Rotate_x(-time_delta * CAMERA_ROTATE_SPEED);
}

void RotateCamraLeft(double time_delta)
{
    view_frame.normal = view_frame.normal.Rotate_y(-time_delta * CAMERA_ROTATE_SPEED);
}

void RotateCameraRight(double time_delta)
{
    view_frame.normal = view_frame.normal.Rotate_y(time_delta * CAMERA_ROTATE_SPEED);
}

void RotateCameraDownwards(double time_delta)
{
    view_frame.normal = view_frame.normal.Rotate_x(time_delta * CAMERA_ROTATE_SPEED);
}

void KeyInput(std::chrono::high_resolution_clock::time_point& time_last)
{
    std::chrono::duration<double> time_delta = std::chrono::high_resolution_clock::now() - time_last;

    // WASD keys
    if(W_KEY_PRESSED)
        MoveCameraForward(time_delta.count());
    if(A_KEY_PRESSED)
        MoveCameraLeft(time_delta.count());
    if(S_KEY_PRESSED)
        MoveCameraBackward(time_delta.count());
    if(D_KEY_PRESSED)
        MoveCameraRight(time_delta.count());

    // Arrow keys
    if(UP_KEY_PRESSED)
        RotateCameraUpwards(time_delta.count());
    if(LEFT_KEY_PRESSED)
        RotateCamraLeft(time_delta.count());
    if(RIGHT_KEY_PRESSED)
        RotateCameraRight(time_delta.count());
    if(DOWN_KEY_PRESSED)
        RotateCameraDownwards(time_delta.count());
}

void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods)
{
    // Escape key
    if (key == GLFW_KEY_ESCAPE && action == GLFW_PRESS)
        exit(0);

    // WASD keys
    else if(key == GLFW_KEY_W && action == GLFW_PRESS)
        W_KEY_PRESSED = true;
    else if(key == GLFW_KEY_W && action == GLFW_RELEASE)
        W_KEY_PRESSED = false;

    else if(key == GLFW_KEY_A && action == GLFW_PRESS)
        A_KEY_PRESSED = true;
    else if(key == GLFW_KEY_A && action == GLFW_RELEASE)
        A_KEY_PRESSED = false;

    else if(key == GLFW_KEY_S && action == GLFW_PRESS)
        S_KEY_PRESSED = true;
    else if(key == GLFW_KEY_S && action == GLFW_RELEASE)
        S_KEY_PRESSED = false;

    else if(key == GLFW_KEY_D && action == GLFW_PRESS)
        D_KEY_PRESSED = true;
    else if(key == GLFW_KEY_D && action == GLFW_RELEASE)
        D_KEY_PRESSED = false;

    
    // Arrow keys
    else if(key == GLFW_KEY_UP && action == GLFW_PRESS)
        UP_KEY_PRESSED = true;
    else if(key == GLFW_KEY_UP && action == GLFW_RELEASE)
        UP_KEY_PRESSED = false;

    else if(key == GLFW_KEY_LEFT && action == GLFW_PRESS)
        LEFT_KEY_PRESSED = true;
    else if(key == GLFW_KEY_LEFT && action == GLFW_RELEASE)
        LEFT_KEY_PRESSED = false;

    else if(key == GLFW_KEY_RIGHT && action == GLFW_PRESS)
        RIGHT_KEY_PRESSED = true;
    else if(key == GLFW_KEY_RIGHT && action == GLFW_RELEASE)
        RIGHT_KEY_PRESSED = false;

    else if(key == GLFW_KEY_DOWN && action == GLFW_PRESS)
        DOWN_KEY_PRESSED = true;
    else if(key == GLFW_KEY_DOWN && action == GLFW_RELEASE)
        DOWN_KEY_PRESSED = false;
}

/**************************************************************************************************************/


/*
    Return the index of the entity in "entities" that the given ray intersects first, -1 if there are no intersections.
    If the optional last argument is provided, the given index into "entities" will be ignored during the search.
*/
int NearestIntersection(Ray& ray, Vec3<double>& intersection, double& distance, int entity_index_ignore = -1)
{
    double temp_distance;
    Vec3<double> temp_intersection;
    int nearest_entity_index = -1;
    distance = std::numeric_limits<double>::max(); // Start with an impossibly large distance to compare with

    // Go through all entities in the scene
    for(int i = 0; i < entities.size(); ++i)
    {
        // If this entity is not flagged to be ignored, and the ray has intersected it
        if(i != entity_index_ignore && entities[i]->RayIntersect(ray, temp_intersection))
        {
            // Check if this intersection is closer than one we have already found
            temp_distance = ray.origin.Distance(temp_intersection);

            // If so, mark it as such
            if(temp_distance < distance)
            {
                intersection = temp_intersection;
                distance = temp_distance;
                nearest_entity_index = i;
            }
        }
    }

    return nearest_entity_index;
}


/*
    Scale a Color given a lighting adjustment value.
    Extra checking must be done, as Colors are stored as unsigned chars. Color values are temporarily stored as ints and then converted back to unsigned chars
    to avoid overflow/underflow.
*/
void AdjustColorLighting(Color& c, double adjustment)
{
    // Adjust Color value, ensuring there is no overflow/underflow
    auto adjust_and_check = [&adjustment](GLubyte& value){
        int temp = value;
        temp *= adjustment;
        if(temp < 0)
            temp = 0;
        else if(temp > COLOR_VALUE_MAX)
            temp = COLOR_VALUE_MAX;
        
        value = temp;
    };

    adjust_and_check(c.b);
    adjust_and_check(c.g);
    adjust_and_check(c.r);
}

/*
    Find the lighting value at a point.
    The index of the lit entity is provided so that it can be ignored when checking paths between it and the lights. Since the point of incidence
    is itself on the lit_entity, it could possibly be considered to intersect the latter. We avoid this by ignoring the lit_entity_index.
*/
double RayLight(Vec3<double> incidence, Ray normal_at_incident, int lit_entity_index)
{
    double light_additive = 0.0;

    Ray light_ray;
    Vec3<double> intersection;
    int entity_index;
    double distance;

    light_ray.origin = incidence;

    // Check each light to see if there is an uninterrupted path between it and the entity
    for(int i = 0; i < lights.size(); ++i)
    {
        light_ray.direction = lights[i].origin - incidence;
        light_ray.direction = light_ray.direction.Normalize();


        // If the light_ray intersects some entity
        if((entity_index = NearestIntersection(light_ray, intersection, distance, lit_entity_index)) != -1)
            continue;

        light_additive += lights[i].intensity * std::max(0.0, (light_ray.direction.DotProduct(normal_at_incident.direction)));

    }

    return light_additive;
}

/*
    Find the Color value associated with the given ray.
*/
Color RayColor(Ray& r, int depth, double reflected_light_factor = 0.0)
{
    Color c = BACKGROUND_COLOR;
    double distance = std::numeric_limits<double>::max(); // Start with an impossibly large distance to compare with
    Vec3<double> intersection;
    int entity_index;

    // If an entity was intersected by the ray
    if((entity_index = NearestIntersection(r, intersection, distance)) != -1)
    {
        Entity* intersected_entity = entities[entity_index].get();
        if(intersected_entity->surface == Surface_Type::OPAQUE)
        {
            // Get the color of the entity that we intersected
            c = intersected_entity->color;

            // Find the lighting factor at the intersection point
            Ray normal_at_intersection = intersected_entity->NormalAtPoint(intersection);
            double light_additive = RayLight(intersection, normal_at_intersection, entity_index);
            AdjustColorLighting(c, light_additive + reflected_light_factor);

        }
        else if(intersected_entity->surface == Surface_Type::REFLECTIVE)
        {
            if(depth == REFLECTIVE_DEPTH_MAX)
            {
                c = intersected_entity->color;
            }
            else
            {
                // Calculate the reflection ray
                Ray reflected_ray = intersected_entity->ReflectedRay(r, intersection);

                // Calculate the lighting factor at the intersection point.
                Ray normal_at_intersection = intersected_entity->NormalAtPoint(intersection);
                double light_factor = RayLight(intersection, normal_at_intersection, entity_index);

                // Call RayColor recursively to find the color associated with the ray that was reflected.
                Color reflection_color = RayColor(reflected_ray, depth+1, intersected_entity->reflectivity * light_factor);
                c = (1.0 - intersected_entity->reflectivity) * intersected_entity->color + intersected_entity->reflectivity * reflection_color;
                AdjustColorLighting(c, light_factor);
            }
        }
    }

    return c;
}

/*
    Struct that contains the information that is passed on to each thread.
*/
struct thread_info{
    unsigned int id;                        // Id of the thread
    Vec3<double> bottom_left_pixel_center;  // Coord of the bottom left pixel of the viewframe
    Vec3<double> q_x, q_y;                  // Amounts to shift in the x and y directions to get to the next pixel
};

/*
    Entry point for threads.
    Each thread will continuously calculate the colors associated with each ray that is shot through a viewframe pixel.
*/
void RayTrace_thread(void *info)
{
    thread_info ti = *(thread_info*)info;

    unsigned int stride = thread_pool.size();

    unsigned int x, y;
    unsigned int aligned_index;
    Color c;

    Vec3<double> point;
    Ray ray;

    // Calculate the color for each pixel until there are no more, incrementing by "stride" amount of pixels each loop.
    for(int i = ti.id; i < PIXEL_DATA_SIZE / COLOR_DEPTH; i += stride)
    {
        // Convert i, which represents the "flat index" of the pixel array into 2-D coordinates.
        x = i % VIEW_FRAME_WIDTH;
        y = i / VIEW_FRAME_WIDTH;
        aligned_index = i * COLOR_DEPTH; // Each pixel is represented by 3 contiguous indexes in the pixel_data array (one index for each color: BGR)

        point = camera + ti.bottom_left_pixel_center + (ti.q_x * x) + (ti.q_y * y);

        ray.origin = camera;
        ray.direction = point - camera;
        ray.direction = ray.direction.Normalize();

        c = RayColor(ray, 0);

        pixel_data[aligned_index] = c.b;
        pixel_data[aligned_index + 1] = c.g;
        pixel_data[aligned_index + 2] = c.r;
    }

    free(info);
}

/*
    Starting point of the actual ray tracing.
*/
void RayTrace()
{   
    Vec3<double> viewport_center = camera + (view_frame.offset * view_frame.normal.Normalize());    // Center of our viewpower, given the current viewframe
    Vec3<double> target_vec = viewport_center - camera; // Direction of the ray pointing from the camera to the viewport center
    Vec3<double> b = roll_component.CrossProduct(target_vec);   // Vector represeting the positive x direction when facing the viewport center from the camera
    Vec3<double> target_vec_normal = target_vec.Normalize();
    Vec3<double> b_normal = b.Normalize();
    Vec3<double> v_normal = target_vec_normal.CrossProduct(b_normal);   // Vector represeting the positive y direction when facing the viewport center from the camera

    double g_x = view_frame.offset * std::tan(FOV / 2);                 // Half the size of the viewport's width
    double g_y = g_x * VIEW_FRAME_HEIGHT / (double)VIEW_FRAME_WIDTH;    // Half the size of the viewport's height

    Vec3<double> q_x = (2 * g_x / (double)(VIEW_FRAME_WIDTH - 1)) * b_normal;   // Pixel shifting vector along the width
    Vec3<double> q_y = (2 * g_y / (double)(VIEW_FRAME_HEIGHT - 1)) * v_normal;  // Pixel shifting vector along the height

    Vec3<double> bottom_left_pixel_center = (target_vec_normal * view_frame.offset) - (g_x * b_normal) - (g_y * v_normal);

    // Spawn threads and assign their tasks
    thread_info *ti;

    for(unsigned int i = 0; i < thread_pool.size(); ++i)
    {
        ti = new thread_info{i, bottom_left_pixel_center, q_x, q_y};
        thread_pool[i] = std::thread(RayTrace_thread, (void*)ti);
    }

    for(std::thread &t : thread_pool)
    {
        t.join();
    }

    // The pixel_data is completed for this frame, so give it to OpenGL to draw
    glDrawPixels(VIEW_FRAME_WIDTH, VIEW_FRAME_HEIGHT, GL_BGR_EXT, GL_UNSIGNED_BYTE, pixel_data);
}

void UpdateFPS(GLFWwindow *window)
{
    static double last_time = glfwGetTime();
    static int frame_count = 0;

    double current_time = glfwGetTime();
    double delta = current_time - last_time;

    frame_count++;
    if(delta >= 1.0)
    {
        double fps = static_cast<double>(frame_count) / delta;

        std::stringstream ss;
        ss << WINDOW_TITLE << " - " << fps << " fps";

        glfwSetWindowTitle(window, ss.str().c_str());

        frame_count = 0;
        last_time = current_time;
    }
}

void InitPlainFloor()
{
    entities.push_back(
        std::unique_ptr<Entity>( new
            Rectangle_Entity{   {-100.0, -10.0, 1000.0},    // Point A
                                {100.0, -10.0, 1000.0},     // Point B
                                {-100.0, -10.0, -100.0},    // Point C
                                
                                {50, 50, 100},          // Color
                                Surface_Type::OPAQUE,   // Surface type
                                0.0})                   // Reflectivity
        );
}

void InitSceneEntities()
{
    
    // Add spheres
    entities.push_back( 
        std::unique_ptr<Entity>( new

        Sphere_Entity{  {{-4.0, 2.0, 60.0}, {6.0}}, // sphere.origin, sphere.radius
                         {125, 125, 125},           // color
                          Surface_Type::REFLECTIVE, // surface type
                           0.5}                     // reflectivity
    ));

    entities.push_back( 
        std::unique_ptr<Entity>( new

        Sphere_Entity{  {{0.0, 12.0, 55.0}, {4.0}}, // sphere.origin, sphere.radius
                         {125, 125, 0},             // color
                          Surface_Type::REFLECTIVE, // surface type
                           0.8}                     // reflectivity
    ));

    entities.push_back(
        std::unique_ptr<Entity>( new

        Sphere_Entity{  {{5.0, 0.0, 50.0}, {4.0}},  // sphere.origin, sphere.radius
                         {125, 50, 0},              // color
                          Surface_Type::OPAQUE,     // surface type
                          0.0}                      // reflectivity
    ));
    
    // Add floor
    InitPlainFloor();
}

int main()
{
    std::chrono::high_resolution_clock::time_point time_last;
    GLFWwindow* window;

    /* Initialize the library */
    if (!glfwInit())
        return -1;

    /* Create a windowed mode window and its OpenGL context */
    window = glfwCreateWindow(VIEW_FRAME_WIDTH, VIEW_FRAME_HEIGHT, WINDOW_TITLE.data(), NULL, NULL);
    if (!window)
    {
        glfwTerminate();
        return -1;
    }

    /* Make the window's context current */
    glfwMakeContextCurrent(window);

    // Enable keyboard input
    glfwSetKeyCallback(window, key_callback);

    // Specify the alignment requirements for the start of each pixel row in memory
    glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
    glPixelStorei(GL_PACK_ALIGNMENT, 1);

    // Allocate thread pool
    thread_pool = std::vector<std::thread>(std::thread::hardware_concurrency() - 1); // Subtract one to account for main thread, which will manage spawned threads

    InitSceneEntities();

    time_last = std::chrono::high_resolution_clock::now();

    /* Loop until the user closes the window */
    while (!glfwWindowShouldClose(window))
    {
        /* Render here */
        glClearColor(0.0, 0.0, 0.0, 1);
        glClear(GL_COLOR_BUFFER_BIT);

        RayTrace();

        /* Swap front and back buffers */
        glfwSwapBuffers(window);

        /* Poll for and process events */
        glfwPollEvents();
        KeyInput(time_last);
        time_last = std::chrono::high_resolution_clock::now();

        UpdateFPS(window);
    }

    glfwTerminate();
    return 0;
}