# AtomMap

AtomMap Version 1.1 :
---------------------

AtomMap maps atoms in chemical reactions.

Installation
------------

AtomMap requires OCaml compiler (version 4.02.3 or newer) together with ``xml-light`` library
(package ``libxml-light-ocaml-dev`` for Ubuntu) and ``xlib`` library.

In order to install ``xlib`` library clone:

    git clone https://github.com/w-jaworski/xlib

and install it according to README provided.

Compile AtomMap:

    make

mkdir results (dodać do makefile)

The following executable files will be created:
``map_reaction_server`` establish a server 

``map_reaction_server`` is a socket server that awaits for queries on port 2727.
Queries are reactions in SMILES format to be mapped.
Matching results are provded in XML format.

Credits
-------

Copyright (C) 2015-2017 Wojciech Jaworski <wjaworski atSPAMfree mimuw dot edu dot pl>
Copyright (C) 2015-2017 Institute of Informatics, University of Warsaw 

Copyright (C) 2015-2017 Sara Szymkuc <saraszymkuc atSPAMfree gmail dot com>
Copyright (C) 2015-2017 Barbara Mikulak <basia dot mikulak atSPAMfree gmail dot com>
Copyright (C) 2015-2017 Institute of Organic Chemistry, Polish Academy of Sciences

The economic rights to the work belongs in 50% to University of Warsaw and
in 50% to Institute of Organic Chemistry, Polish Academy of Sciences
