CXX=g++

CFlags=-c -O3 -MMD -std=c++17
Sources=$(wildcard src/*.cpp src/examples/*.cpp)
IncludeDir=-I./include
AllObjects=$(patsubst src/%, obj/%, $(Sources:.cpp=.o))
Executables=main
Objects=$(filter-out $(addprefix obj/,$(Executables:=.o)),$(AllObjects))

all: $(Sources) $(Executables)

$(Executables): $(AllObjects)
	@mkdir -p data obj obj/examples
	$(CXX) $(Objects) $(addprefix obj/,$@.o) $(GSLFlags) -o $@

obj/%.o: src/%.cpp
	@mkdir -p $(@D)
	$(CXX) $(CFlags) $(IncludeDir) $< -o $@

-include $(AllObjects:.o=.d)

clean:
	rm -rf obj/*.o obj/*.d obj/examples/*.o obj/examples/*.d $(Executables)
