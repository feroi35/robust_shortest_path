#include "parser.h"


Instance::Instance(IloEnv env, char filename[]) {
    // I'm not sure why the instance depends of the environment
    // and I think it make more sense not to do so
    // But that's what the professor did in the TP
    name = filename;
    p = IloNumArray(env);
    ph = IloNumArray(env);
    mat = std::vector<Arc>();
    
    char readChar;
    int readInt;
    std::ifstream file(filename);

    if (!file) {
        std::cerr << "Error when opening file" << std::endl;
        exit(1);
    }

    file >> readChar >> readChar >> n;
    file >> readChar >> readChar >> s;
    file >> readChar >> readChar >> t;
    file >> readChar >> readChar >> S;
    file >> readChar >> readChar >> readChar >> d1;
    file >> readChar >> readChar >> readChar >> d2;
    file >> readChar >> readChar;
    for (unsigned int k=0; k<n; k++) {
        file >> readChar >> readInt;
        p.add(readInt);
    }
    // ] ph =
    file >> readChar >> readChar >> readChar >> readChar;
    for (unsigned int k=0; k<n; k++) {
        file >> readChar >> readInt;
        ph.add(readInt);
    }
    // ] mat = [
    file >> readChar >> readChar >> readChar >> readChar >> readChar >> readChar;

    std::vector<std::vector<double>> d_(n, std::vector<double>(n, undefinedValue));
    std::vector<std::vector<double>> D_(n, std::vector<double>(n, undefinedValue));
    // débattable, c'est pour avoir une première version qui tourne

    while (readChar != ']') {
        Arc v;
        file >> v.i;
        file >> v.j;
        file >> v.d;
        file >> v.D;
        file >> readChar; // either ';' or ']'
        mat.push_back(v);
        d_[v.i-1][v.j-1] = v.d;
        D_[v.i-1][v.j-1] = v.D;
    }
    file.close();
    d = d_;
    D = D_;
}


void Instance::display() const {
    std::cout << "n = " << n << std::endl;
    std::cout << "s = " << s << std::endl;
    std::cout << "t = " << t << std::endl;
    std::cout << "S = " << S << std::endl;
    std::cout << "d1 = " << d1 << std::endl;
    std::cout << "d2 = " << d2 << std::endl;
    std::cout << "p = [";
    for (unsigned int i=0; i<n; i++) {
        std::cout << p[i] << " ";
    }
    std::cout << "]" << std::endl;
    std::cout << "ph = [";
    for (unsigned int i=0; i<n; i++) {
        std::cout << ph[i] << " ";
    }
    std::cout << "]" << std::endl;

    std::cout << "mat = [" << std::endl;
    for (unsigned int i=0; i<mat.size(); i++) {
        std::cout << "(" << mat[i].i << ", " << mat[i].j << ", " << mat[i].d << ", " 
            << mat[i].D << ") " << std::endl;
    }
    std::cout << "]" << std::endl;
    std::cout << "d = [" << std::endl;
    for (unsigned int i=0; i<d.size(); i++) {
        std::cout << "[";
        for (unsigned int j=0; j<d[i].size(); j++) {
            std::cout << d[i][j] << " ";
        }
        std::cout << "]" << std::endl;
    }
}

double Instance::compute_static_score(std::vector<IloInt> sol) const {
    double static_score = 0.0;
    IloInt current_node = sol[0];
    assert(current_node == s); 
    for (unsigned int i=1; i<sol.size(); i++) {
        IloInt next_node = sol[i];
        static_score += d[current_node][next_node];
        current_node = next_node;
    }
    assert(current_node == t);
    return static_score;
}
