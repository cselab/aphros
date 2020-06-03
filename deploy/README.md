# Deploy

*   Select profile and create `ap.setenv`

     ```
     ./install_setenv PREFIX
     (daint) ./install_setenv PREFIX --profile daint
     ```
*   Set environment

     ```
     . ap.setenv
     ```

*   Install libraries and tools

     ```
     mkdir build
     cd build
     cmake ..
     make -j4
     make install
     ```

     Alternatively, use `ccmake ..` to configure with a dialog.
