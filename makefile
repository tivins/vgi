objs=obj/tuple.o obj/util.o obj/geom.o obj/scene.o obj/trace.o
opts=-Wall -Wextra
bins=vgit

all: $(bins)

vgit: $(objs) obj/main.o
	g++ -Wall -Wextra -o $@ $^

obj/%.o: src/%.cpp
	g++ -std=c++11 $(opts) -c $^ -o $@

clean:
	$(RM) $(objs)

proper: clean
	$(RM) $(bins)