g++ main.cpp -o sfml-window \
                                                                          -I/Users/mario1234/Downloads/SFML-3.0.2/include \
                                                                          -L/Users/mario1234/Downloads/SFML-3.0.2/lib \
                                                                          -lsfml-graphics -lsfml-window -lsfml-system -std=c++17
DYLD_LIBRARY_PATH=/Users/mario1234/Downloads/SFML-3.0.2/lib ./sfml-window