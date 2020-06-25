CC = g++
CFLAGS = -std=c++17 -lglfw3 -lopengl32 -lgdi32 -O2

all: main

main: main.cpp ray_math.cpp
	$(CC) main.cpp ray_math.cpp -o raytracer $(CFLAGS)

clean:
	rm *.exe