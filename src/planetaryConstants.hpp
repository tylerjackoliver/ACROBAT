#include <unordered_map>
#include <vector>

std::unordered_map<std::string, double> PLANETARY_RADII({{"Sun", 1.0},
                                                        {"Earth", 1.0},
                                                        {"Mercury", 1.0},
                                                        {"Jupiter", 1.0}});

std::unordered_map<std::string, double> PLANETARY_GM    ({{"Sun", 1.32712440042e11},
                                                        {"Mercury", 0.022032e6},
                                                        {"Jupiter", 1.266865341960128E+08},
                                                        {"Saturn", 3.793120615901047E+07},
                                                        {"Pluto", 8.693390780381926E+02},
                                                        {"Charon", 1.062509269522026E+02},
                                                        {"Io", 5.959924010272514E+03},
                                                        {"Callisto", 7.179304867611079E+03},
                                                        {"Ganymede", 9.887819980080976E+03},
                                                        {"Triton", 1.428964197232701E+03},
                                                        {"Titan", 8.978137369591670E+03},
                                                        {"Saturn", 3.793120615901047E+07},
                                                        {"Europa", 3.202739815114734E+03},
                                                        {"Venus", 324858.592000},
                                                        {"Moon", 4.9048695e3},
                                                        {"Earth",  3.986004418e5}});