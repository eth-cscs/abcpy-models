// This file is intended to write a bitmap image file corresponding to
// the output produced from the heatedplate.c or heatedplate.f90
// program files.
//
// The BmpImage class was written by Daniel LePage, 2007
// The main method was written by Aaron Bloomfield, 2008

# include "grid_to_bmp.hpp"

int main ( int argc, char** argv ){
    cout << "\n";
    cout << "GRID_TO_BMP:\n";
    cout << "  C++ version\n";

    if ( argc != 3 )
    {
        cerr << "\n";
        cerr << "GRID_TO_BMP\n";
        cerr << "  USAGE: " << argv[0] << " <input-filename> <output-filename>\n" << endl;
        return 1;
    }
    //
    ImageWriter w(argv[1], argv[2]);
    w.write();
}
