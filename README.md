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

The following executable files will be created:
``map_reaction`` maps single reaction
``map_reaction_server`` establishes a server 

``map_reaction`` takes a single command line argument which is a reaction in SMILES format
and generate matching results are provded in XML format on standart output. For example:

    ./map_reaction "ClC(C1=CC=CC=C1)=O.C1=CC=CC=C1>>O=C(C1=CC=CC=C1)C2=CC=CC=C2.[H]Cl"

``map_reaction_server`` is a socket server that awaits for queries on port 2727.
A new process is created for each connection.
Queries are reactions in SMILES format to be mapped.
Matching results are provded in XML format. For example:

    ./map_reaction_server &
    echo "ClC(C1=CC=CC=C1)=O.C1=CC=CC=C1>>O=C(C1=CC=CC=C1)C2=CC=CC=C2.[H]Cl" | netcat localhost 2727

Credits
-------

Copyright (C) 2015-2017 Wojciech Jaworski <wjaworski atSPAMfree mimuw dot edu dot pl>
Copyright (C) 2015-2017 Institute of Informatics, University of Warsaw 

Copyright (C) 2015-2017 Sara Szymkuc <saraszymkuc atSPAMfree gmail dot com>
Copyright (C) 2015-2017 Barbara Mikulak <basia dot mikulak atSPAMfree gmail dot com>
Copyright (C) 2015-2017 Institute of Organic Chemistry, Polish Academy of Sciences

The economic rights to the work belongs in 50% to University of Warsaw and
in 50% to Institute of Organic Chemistry, Polish Academy of Sciences
