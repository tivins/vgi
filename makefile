


all:
	g++ -Wall -Wextra -o a.exe src/main.cpp -I libs/sdl2/include -L libs/sdl2/lib -lSDL2main -lSDL2