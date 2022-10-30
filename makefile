OCAMLC=ocamlc
OCAMLOPT=ocamlopt
OCAMLDEP=ocamldep
#INCLUDES=-I +xml-light -I +gsl -I +xlib -I +zip -I +bz2 -I +zarith -I /usr/local/lib/ocaml/4.01.0/ANSITerminal 
INCLUDES=-I +xml-light -I +gsl -I +xlib -I +zip -I +bz2 -I +zarith
OCAMLFLAGS=$(INCLUDES)
#OCAMLOPTFLAGS=$(INCLUDES) unix.cmxa xml-light.cmxa str.cmxa bigarray.cmxa nums.cmxa zip.cmxa bz2.cmxa xlib.cmxa ANSITerminal.cmxa
OCAMLOPTFLAGS=$(INCLUDES) unix.cmxa xml-light.cmxa str.cmxa bigarray.cmxa zip.cmxa bz2.cmxa zarith.cmxa xlib.cmxa

MACZOWANIE= collection.mli collection.ml pair.mli pair.ml labels.mli labels.ml types.ml hash.mli hash.ml isomorphism.mli isomorphism.ml commonSubstructure.mli commonSubstructure.ml

MACZOWANIE2= $(MACZOWANIE) smiles.ml patterns.ml reactionClasses.ml atomMapping.ml import.ml molecule.ml matchingExec.ml


all:
	# $(OCAMLOPT) -o cr $(OCAMLOPTFLAGS) $(MACZOWANIE2) matching_tests.ml
	$(OCAMLOPT) -o map_reaction $(OCAMLOPTFLAGS) $(MACZOWANIE2) mapReaction.ml
#	$(OCAMLOPT) -o map_reactions $(OCAMLOPTFLAGS) $(MACZOWANIE2) mapReactions.ml
#	$(OCAMLOPT) -o chem.cgi $(OCAMLOPTFLAGS) collection.mli collection.ml pair.mli pair.ml labels.mli labels.ml types.ml smiles.ml import.ml chem_cgi.ml
	$(OCAMLOPT) -o map_reaction_server $(OCAMLOPTFLAGS) $(MACZOWANIE2) server.ml
#	$(OCAMLOPT) -o map_reaction_worker $(OCAMLOPTFLAGS) $(MACZOWANIE2) worker.ml
	mkdir -p results

.SUFFIXES: .mll .mly .ml .mli .cmo .cmi .cmx

.mll.ml:
	ocamllex $<

.mly.mli:
	ocamlyacc $<

.mly.ml:
	ocamlyacc $<

.ml.cmo:
	$(OCAMLC) $(OCAMLFLAGS) -c $<

.mli.cmi:
	$(OCAMLC) $(OCAMLFALGS) -c $<

.ml.cmx:
	$(OCAMLOPT) $(OCAMLOPTFLAGS) -c $<

xlib.cmxa:
	cd xlib; make $@

xlib.cma:
	cd xlib; make $@

clean:
	rm -f *~ *.cm[oix] *.o cr map_reaction map_reactions chem.cgi map_reaction_server map_reaction_worker
