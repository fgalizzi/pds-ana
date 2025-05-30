// Link reference: https://docs.google.com/spreadsheets/d/1N9xcb2VVlzzDcNfBjlj_buhH9LiBTdG8-cnisb-orsI/edit?gid=926497567#gid=926497567 
// Analysis performed with...

#include <vector>

// M1 (20,27), M2 (21,26), M3 (0,2), M4 (1,3)


// --- M1 ------------------------------------------------------------------
int module = 1;
std::vector<int> module_channels = {20, 27};
std::vector<double> v_brs     = {42.47, 41.84};
std::vector<double> err_v_brs = {0.21, 0.24};
std::vector<int> vgains = {1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900,
                           2000, 2100, 2200, 2300, 2400, 2500, 2600, 2700, 2800, 2900,
                           3000};
double bias_dac = 1200;
std::vector<int> runs   = {34023, 34024, 34025, 34026, 34027, 34028, 34029, 34030, 34031, 34032,
                           34033, 34034, 34035, 34036, 34037, 34038, 34039, 34040, 34041, 34042,
                           34058};

double bias_dac = 1187;
std::vector<int> runs   = {33880, 33881, 33882, 33883, 33884, 33885, 33886, 33887, 33888, 33889,
                           33890, 33891, 33892, 33893, 33894, 33895, 33896, 33897, 33898, 33899,
                           34059};

double bias_dac = 1174;
std::vector<int> runs   = {33900, 33901, 33902, 33903, 33904, 33905, 33906, 33907, 33908, 33909,
                           33910, 33911, 33912, 33913, 33914, 33915, 33916, 33917, 33918, 33919,
                           34060};

double bias_dac = 1161;
std::vector<int> runs   = {33920, 33921, 33922, 33923, 33924, 33925, 33926, 33927, 33928, 33929,
                           33930, 33931, 33932, 33933, 33934, 33935, 33936, 33937, 33938, 33939,
                           34061};

double bias_dac = 1148;
std::vector<int> runs   = {33940, 33941, 33942, 33943, 33944, 33945, 33946, 33947, 33948, 33949,
                           33950, 33951, 33952, 33953, 33954, 33955, 33956, 33957, 33958, 33959,
                           34062};

std::vector<double> bias_dacs = {1200, 1187, 1174, 1161, 1148};
std::vector<std::vector<int>> run_batches = {{34023, 34024, 34025, 34026, 34027, 34028, 34029, 34030, 34031, 34032,
                                      34033, 34034, 34035, 34036, 34037, 34038, 34039, 34040, 34041, 34042,
                                      34058},
                                     {33880, 33881, 33882, 33883, 33884, 33885, 33886, 33887, 33888, 33889,
                                      33890, 33891, 33892, 33893, 33894, 33895, 33896, 33897, 33898, 33899,
                                      34059},
                                     {33900, 33901, 33902, 33903, 33904, 33905, 33906, 33907, 33908, 33909,
                                      33910, 33911, 33912, 33913, 33914, 33915, 33916, 33917, 33918, 33919,
                                      34060},
                                     {33920, 33921, 33922, 33923, 33924, 33925, 33926, 33927, 33928, 33929,
                                      33930, 33931, 33932, 33933, 33934, 33935, 33936, 33937, 33938, 33939,
                                      34061},
                                     {33940, 33941, 33942, 33943, 33944, 33945, 33946, 33947, 33948, 33949,
                                      33950, 33951, 33952, 33953, 33954, 33955, 33956, 33957, 33958, 33959,
                                      34062}};

// --- M2 ------------------------------------------------------------------
int module = 2;
std::vector<int> module_channels = {21, 26};
std::vector<double> v_brs     = {42.71, 42.59};
std::vector<double> err_v_brs = {0.05, 0.06};
std::vector<int> vgains = {0,    100,  200,  300,  400,  500,  600,  700,  800,  900,
                           1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900,
                           2000, 2100, 2200, 2300, 2400, 2500, 2600, 2700, 2800, 2900,
                           3000};

double bias_dac = 1200; 
std::vector<int> runs   = {34230, 34231, 34232, 34233, 34234, 34235, 34236, 34237, 34238, 34239,
                           33380, 33381, 33382, 33383, 33384, 33385, 33386, 33387, 33388, 33389,
                           33390, 33392, 33393, 33394, 33395, 33396, 33397, 33398, 33399, 33400,
                           34043};
std::vector<int> runs   = {34230, 34231, 34232, 34233, 34234, 34235, 34236, 34237, 34238, 34239,
                           33963, 33964, 33965, 33966, 33967, 33968, 33969, 33970, 33971, 33972,
                           33973, 33974, 33975, 33976, 33977, 33978, 33979, 33980, 33981, 33982,
                           34043};

double bias_dac = 1187;
std::vector<int> runs   = {34240, 34241, 34242, 34243, 34244, 34245, 34246, 34247, 34248, 34249,
                           33401, 33402, 33403, 33404, 33405, 33406, 33407, 33408, 33409, 33410,
                           33411, 33413, 33414, 33415, 33416, 33417, 33418, 33419, 33420, 33421, 
                           34044};
std::vector<int> runs   = {34240, 34241, 34242, 34243, 34244, 34245, 34246, 34247, 34248, 34249,
                           33640, 33641, 33642, 33643, 33644, 33645, 33646, 33647, 33648, 33649,
                           33650, 33651, 33652, 33653, 33654, 33655, 33656, 33657, 33658, 33659,
                           34044};

double bias_dac = 1174;
std::vector<int> runs   = {34250, 34251, 34252, 34253, 34254, 34255, 34256, 34257, 34258, 34259,
                           33660, 33661, 33662, 33663, 33464, 33665, 33666, 33667, 33668, 33669, 
                           33670, 33671, 33672, 33673, 33674, 33675, 33676, 33677, 33678, 33679,
                           34045};

double bias_dac = 1161;
std::vector<int> runs   = {34260, 34261, 34262, 34263, 34264, 34265, 34266, 34267, 34268, 34269,
                           33680, 33681, 33682, 33683, 33684, 33685, 33686, 33687, 33688, 33689, 
                           33690, 33691, 33692, 33693, 33694, 33695, 33696, 33697, 33698, 33699,
                           34046};

double bias_dac = 1148;
std::vector<int> runs   = {34270, 34271, 34272, 34273, 34274, 34275, 34276, 34277, 34278, 34279,
                           33700, 33701, 33702, 33703, 33704, 33705, 33706, 33707, 33708, 33709, 
                           33710, 33711, 33712, 33713, 33714, 33715, 33716, 33717, 33718, 33719,
                           34047};

std::vector<double> bias_dacs = {1200, 1187, 1174, 1161, 1148};
std::vector<std::vector<int>> run_batches = {{34230, 34231, 34232, 34233, 34234, 34235, 34236, 34237, 34238, 34239,
                                      33380, 33381, 33382, 33383, 33384, 33385, 33386, 33387, 33388, 33389,
                                      33390, 33392, 33393, 33394, 33395, 33396, 33397, 33398, 33399, 33400,
                                      34043},
                                     {34240, 34241, 34242, 34243, 34244, 34245, 34246, 34247, 34248, 34249,
                                      33401, 33402, 33403, 33404, 33405, 33406, 33407, 33408, 33409, 33410,
                                      33411, 33413, 33414, 33415, 33416, 33417, 33418, 33419, 33420, 33421, 
                                      34044},
                                     {34250, 34251, 34252, 34253, 34254, 34255, 34256, 34257, 34258, 34259,
                                      33660, 33661, 33662, 33663, 33464, 33665, 33666, 33667, 33668, 33669, 
                                      33670, 33671, 33672, 33673, 33674, 33675, 33676, 33677, 33678, 33679,
                                      34045},
                                     {34260, 34261, 34262, 34263, 34264, 34265, 34266, 34267, 34268, 34269,
                                      33680, 33681, 33682, 33683, 33684, 33685, 33686, 33687, 33688, 33689, 
                                      33690, 33691, 33692, 33693, 33694, 33695, 33696, 33697, 33698, 33699,
                                      34046},
                                     {34270, 34271, 34272, 34273, 34274, 34275, 34276, 34277, 34278, 34279,
                                      33700, 33701, 33702, 33703, 33704, 33705, 33706, 33707, 33708, 33709, 
                                      33710, 33711, 33712, 33713, 33714, 33715, 33716, 33717, 33718, 33719,
                                      34047}};

// --- M3 ------------------------------------------------------------------
int module = 3;
std::vector<int> module_channels = {0, 2};
std::vector<double> v_brs     = {27.45, 27.40};
std::vector<double> err_v_brs = {0.16, 0.16};
std::vector<int> vgains = {0,    100,  200,  300,  400,  500,  600,  700,  800,  900,
                           1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900,
                           2000, 2100, 2200, 2300, 2400, 2500, 2600, 2700, 2800, 2900,
                           3000};


double bias_dac = 806;
std::vector<int> runs   = {34180, 34181, 34182, 34183, 34184, 34185, 34186, 34187, 34188, 34189,
                           33983, 33984, 33985, 33986, 33987, 33988, 33989, 33990, 33991, 33992, 
                           33993, 33994, 33995, 33996, 33997, 33998, 33999, 34000, 34001, 34002,
                           34048};

double bias_dac = 793;
std::vector<int> runs   = {34190, 34191, 34192, 34193, 34194, 34195, 34196, 34197, 34198, 34199,
                           33720, 33721, 33722, 33723, 33724, 33725, 33726, 33727, 33728, 33729, 
                           33730, 33731, 33732, 33733, 33734, 33735, 33736, 33737, 33738, 33739,
                           34049};

double bias_dac = 780;
std::vector<int> runs   = {34200, 34201, 34202, 34203, 34204, 34205, 34206, 34207, 34208, 34209,
                           33740, 33741, 33742, 33743, 33744, 33745, 33746, 33747, 33748, 33749, 
                           33750, 33751, 33752, 33753, 33754, 33755, 33756, 33757, 33758, 33759,
                           34050};

double bias_dac = 767;
std::vector<int> runs   = {34210, 34211, 34212, 34213, 34214, 34215, 34216, 34217, 34218, 34219,
                           33760, 33761, 33762, 33763, 33764, 33765, 33766, 33767, 33768, 33769, 
                           33770, 33771, 33772, 33773, 33774, 33775, 33776, 33777, 33778, 33779,
                           34051};

double bias_dac = 754;
std::vector<int> runs   = {34220, 34221, 34222, 34223, 34224, 34225, 34226, 34227, 34228, 34229,
                           33780, 33781, 33782, 33783, 33784, 33785, 33786, 33787, 33788, 33789, 
                           33790, 33791, 33792, 33793, 33794, 33795, 33796, 33797, 33798, 33799,
                           34052};

std::vector<double> bias_dacs = {806, 793, 780, 767, 754};
std::vector<std::vector<int>> run_batches = {{34180, 34181, 34182, 34183, 34184, 34185, 34186, 34187, 34188, 34189,
                                      33983, 33984, 33985, 33986, 33987, 33988, 33989, 33990, 33991, 33992, 
                                      33993, 33994, 33995, 33996, 33997, 33998, 33999, 34000, 34001, 34002,
                                      34048},
                                     {34190, 34191, 34192, 34193, 34194, 34195, 34196, 34197, 34198, 34199,
                                      33720, 33721, 33722, 33723, 33724, 33725, 33726, 33727, 33728, 33729, 
                                      33730, 33731, 33732, 33733, 33734, 33735, 33736, 33737, 33738, 33739,
                                      34049},
                                     {34200, 34201, 34202, 34203, 34204, 34205, 34206, 34207, 34208, 34209,
                                      33740, 33741, 33742, 33743, 33744, 33745, 33746, 33747, 33748, 33749, 
                                      33750, 33751, 33752, 33753, 33754, 33755, 33756, 33757, 33758, 33759,
                                      34050},
                                     {34210, 34211, 34212, 34213, 34214, 34215, 34216, 34217, 34218, 34219,
                                      33760, 33761, 33762, 33763, 33764, 33765, 33766, 33767, 33768, 33769, 
                                      33770, 33771, 33772, 33773, 33774, 33775, 33776, 33777, 33778, 33779,
                                      34051},
                                     {34220, 34221, 34222, 34223, 34224, 34225, 34226, 34227, 34228, 34229,
                                      33780, 33781, 33782, 33783, 33784, 33785, 33786, 33787, 33788, 33789, 
                                      33790, 33791, 33792, 33793, 33794, 33795, 33796, 33797, 33798, 33799,
                                      34052}};

// --- M4 ------------------------------------------------------------------
int module = 4;
std::vector<int> module_channels = {1, 3};
std::vector<double> v_brs     = {27.42, 27.34};
std::vector<double> err_v_brs = {0.16, 0.17};
std::vector<int> vgains = {0,    100,  200,  300,  400,  500,  600,  700,  800,  900,
                           1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900,
                           2000, 2100, 2200, 2300, 2400, 2500, 2600, 2700, 2800, 2900,
                           3000};
double bias_dac 806;
std::vector<int> runs   = {34130, 34131, 34132, 34133, 34134, 34135, 34136, 34137, 34138, 34139,
                           34003, 34004, 34005, 34006, 34007, 34008, 34009, 34010, 34011, 34012, 
                           34013, 34014, 34015, 34016, 34017, 34018, 34019, 34020, 34021, 34022,
                           34053};

double bias_dac = 793;
std::vector<int> runs   = {34140, 34141, 34142, 34143, 34144, 34145, 34146, 34147, 34148, 34149,
                           33800, 33801, 33802, 33803, 33804, 33805, 33806, 33807, 33808, 33809, 
                           33810, 33811, 33812, 33813, 33814, 33815, 33816, 33817, 33818, 33819,
                           34054};

double bias_dac = 780;
std::vector<int> runs   = {34150, 34151, 34152, 34153, 34154, 34155, 34156, 34157, 34158, 34159,
                           33820, 33821, 33822, 33823, 33824, 33825, 33826, 33827, 33828, 33829, 
                           33830, 33831, 33832, 33833, 33834, 33835, 33836, 33837, 33838, 33839,
                           34055};

double bias_dac = 767;
std::vector<int> runs   = {34160, 34161, 34162, 34163, 34164, 34165, 34166, 34167, 34168, 34169,
                           33840, 33841, 33842, 33843, 33844, 33845, 33846, 33847, 33848, 33849, 
                           33850, 33851, 33852, 33853, 33854, 33855, 33856, 33857, 33858, 33859,
                           34056};

double bias_dac = 754;
std::vector<int> runs   = {34170, 34171, 34172, 34173, 34174, 34175, 34176, 34177, 34178, 34179,
                           33860, 33861, 33862, 33863, 33864, 33865, 33866, 33867, 33868, 33869, 
                           33870, 33871, 33872, 33873, 33874, 33875, 33876, 33877, 33878, 33879,
                           34057};

std::vector<double> bias_dacs = {806, 793, 780, 767, 754};
std::vector<std::vector<int>> run_batches = {{34130, 34131, 34132, 34133, 34134, 34135, 34136, 34137, 34138, 34139,
                                      34003, 34004, 34005, 34006, 34007, 34008, 34009, 34010, 34011, 34012, 
                                      34013, 34014, 34015, 34016, 34017, 34018, 34019, 34020, 34021, 34022,
                                      34053},
                                     {34140, 34141, 34142, 34143, 34144, 34145, 34146, 34147, 34148, 34149,
                                      33800, 33801, 33802, 33803, 33804, 33805, 33806, 33807, 33808, 33809, 
                                      33810, 33811, 33812, 33813, 33814, 33815, 33816, 33817, 33818, 33819,
                                      34054},
                                     {34150, 34151, 34152, 34153, 34154, 34155, 34156, 34157, 34158, 34159,
                                      33820, 33821, 33822, 33823, 33824, 33825, 33826, 33827, 33828, 33829, 
                                      33830, 33831, 33832, 33833, 33834, 33835, 33836, 33837, 33838, 33839,
                                      34055},
                                     {34160, 34161, 34162, 34163, 34164, 34165, 34166, 34167, 34168, 34169,
                                      33840, 33841, 33842, 33843, 33844, 33845, 33846, 33847, 33848, 33849, 
                                      33850, 33851, 33852, 33853, 33854, 33855, 33856, 33857, 33858, 33859,
                                      34056},
                                     {34170, 34171, 34172, 34173, 34174, 34175, 34176, 34177, 34178, 34179,
                                      33860, 33861, 33862, 33863, 33864, 33865, 33866, 33867, 33868, 33869,
                                      33870, 33871, 33872, 33873, 33874, 33875, 33876, 33877, 33878, 33879, 
                                      34057}};
