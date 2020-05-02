{ pkgs ? import <nixpkgs> {} }:
with pkgs.python36Packages; buildPythonApplication {
  name = "bif";
  propagatedBuildInputs = [ numpy ];
}
