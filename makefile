
objs=obj/tuple.o obj/geom.o obj/scene.o obj/trace.o
opts=-Wall -Wextra

all: vgit

vgit: $(objs) obj/main.o
	g++ -Wall -Wextra -o $@ $^

obj/%.o: src/%.cpp
	g++ -std=c++11 $(opts) -c $^ -o $@