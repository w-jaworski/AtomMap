(*
 *  AtomMap: maps atoms in chemical reactions
 *  Copyright (C) 2015-2017 Wojciech Jaworski <wjaworski atSPAMfree mimuw dot edu dot pl> 
 *  Copyright (C) 2015-2017 Institute of Informatics, University of Warsaw                

 *  Copyright (C) 2015-2017 Sara Szymkuc <saraszymkuc atSPAMfree gmail dot com>          
 *  Copyright (C) 2015-2017 Barbara Mikulak <basia dot mikulak atSPAMfree gmail dot com> 
 *  Copyright (C) 2015-2017 Institute of Organic Chemistry Polish Academy of Sciences     

 *  This library is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This library is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *)

open Types
open Xstd
open Printf

(*let test_reaction_stoi filename reactions =
  File.file_out filename (fun file ->
    Xlist.iter reactions (fun r ->
(*       Printf.printf "%s %s\n%!" r.rxn_id r.reaction_smile; *)
      try
        let map = Smiles.calculate_stoi_reaction r in
        let l = List.sort compare (StringMap.fold map [] (fun l name (v,w) ->
          if v = w then l else (name,v,w) :: l)) in
        let l = Xlist.map l (fun (name,v,w) -> Printf.sprintf "%s: %d>>%d" name v w) in
        Printf.fprintf file "%s\t%s\n" r.reaction_smile (String.concat "  " l)
      with e -> Printf.fprintf file "%s\t%s\n" r.reaction_smile (Printexc.to_string e)))*)

(* Sprawdzenie stoimetrii *)
(* let _ = test_reaction_stoi "results/Mapping_algorithm_train_set1_stoi_test.dat" (Import.load_paths ())  *)
(* let _ = test_reaction_stoi "results/OrganicSynthesesReactions_stoi_test.dat" (Import.load_organic_syntheses ()) *)
(* let _ = test_reaction_stoi "results/200_OrganicSyntheses_Apr7_2015_stoi_test.dat" (Import.load_organic_syntheses2 ())  *)

let rec set_invisible_hydrogens_rec visible rev = function
    Atom(p,l) ->
       let rev = if p.name = "H" && not (IntSet.mem visible p.id) then p.id :: rev else rev in
       Xlist.fold l rev (fun rev (_,a) -> set_invisible_hydrogens_rec visible rev a)
  | Link _ -> rev

let set_invisible_hydrogens visible rev = function
    Atom({name = "H"},[]) -> rev
  | Atom({name = "H"},(_,Atom({name = "H"},_)) :: _) -> rev
  | a -> set_invisible_hydrogens_rec visible rev a

let make_id path id =
  if path = "" then id else path ^ "" ^ id

let print_reactions prefix filename reactions =
  File.file_out filename (fun file ->
    Printf.fprintf file "%s" (Smiles.a0poster_header "a2");
    let _,names = Xlist.fold reactions (1,[]) (fun (n,names) re ->
      Printf.fprintf file "%s {%s}\\\\\n" (make_id re.path re.rxn_id) (Smiles.escape_string re.reaction_smile);
      try
      let r,messages = Smiles.prepare_record_for_matching re [] in
      Xlist.iter messages (fun msg -> Printf.fprintf file "{%s}\\\\\n" (Smiles.escape_string msg));
      let rl = [r](*CommonSubstructure.disambiguate_aroma rl*) in
      Xlist.fold rl (n,names) (fun (n,names) r ->
        let l = Xlist.fold r.reactants [] (fun l m -> set_invisible_hydrogens IntSet.empty l m.tree) in
        let l = Xlist.fold r.products l (fun l m -> set_invisible_hydrogens IntSet.empty l m.tree) in
        let labels = Labels.add_invisible_list r.empty_labels l in
        let name = prefix ^ string_of_int n in
(*        File.file_out ("results/images/" ^ name ^ ".xml") (fun xml_file ->
          Printf.fprintf xml_file "%s\n" (Xml.to_string_fmt (Smiles.reaction_to_xml r [] labels IntMap.empty IntMap.empty)));*)
(*         Printf.printf "%d %d\n" (IntSet.size r.reactant_ids) (IntSet.size r.product_ids); *)
        File.file_out ("results/images/" ^ name ^ "r.gv") (fun gv_file ->
          Smiles.reaction_to_gv gv_file 1 r.ids_reactants r [] labels IntMap.empty IntMap.empty);
        File.file_out ("results/images/" ^ name ^ "p.gv") (fun gv_file ->
          Smiles.reaction_to_gv gv_file r.ids_reactants r.ids_products r [] labels IntMap.empty IntMap.empty);
        Printf.fprintf file "\\includegraphics[scale=0.3]{images/%sr.png}$\\to$\n" name;
        Printf.fprintf file "\\includegraphics[scale=0.3]{images/%sp.png}\\\\\n" name;
        n+1,(name ^ "p"(* ^ ".xml"*)) :: (name ^ "r") ::  names)
      with e -> Printf.fprintf file "{%s}\\\\\n" (Smiles.escape_string (Printexc.to_string e));n,names) in
    Sys.chdir "results/images";
(*     ignore (Sys.command ("phantomjs mapping-dl.js " ^ String.concat " " names)); *)
(*     Xlist.iter names (fun name -> ignore (Sys.command ("neato -Tpng " ^ name ^ ".gv -o " ^ name ^ ".png"))); *)
    Xlist.iter names (fun name -> ignore (Sys.command ("neato -Tpng " ^ name ^ ".gv -o " ^ name ^ ".png")));
    Sys.chdir "../..";
    Printf.fprintf file "%s" Smiles.trailer)




let get_bond g x y =
  let p,l = g.(x) in
  let l,n = Xlist.fold l ([],0) (fun (l,n) (b,p) ->
    if p.id = y then l,Smiles.int_of_bond b else (b,p) :: l, n) in
  g.(x) <- p, List.rev l;
  n

let set_bond g x y n =
  if n = 0 then () else
  let p,l = g.(x) in
  g.(x) <- p, (Smiles.bond_of_int n, fst g.(y)) :: l

let add_bonds g bonds =
  let g = Array.copy g in
  Xlist.iter bonds (fun (x,y,n) ->
    let n = get_bond g x y + n in
    let _ = get_bond g y x in
    set_bond g x y n;
    set_bond g y x n);
  g

let match_bond = function
    Triple,Triple -> true
  | Triple,_ -> false
  | _,Triple -> false
  | Double,Double -> true
  | Double,_ -> false
  | _,Double -> false
  | Aromatic,Aromatic -> true
  | Aromatic,_ -> false
  | _,Aromatic -> false
  | _ -> true

(*let rec find_sigmatropic_kernel graph i0 i found rev = function
    [] -> (if i0 < i then List.rev (i :: rev) else i :: rev) :: found
  | b :: bonds -> Xlist.fold (snd graph.(i)) found (fun found (b2,p) ->
       if match_bond (b,b2) && not (Xlist.mem rev p.id) then find_sigmatropic_kernel graph i0 p.id found (i :: rev) bonds else found)*)

let rec find_path graph i0 i found rev = function
    [] -> (if i0 < i then List.rev (i :: rev) else i :: rev) :: found
  | b :: bonds -> Xlist.fold (snd graph.(i)) found (fun found (b2,p) ->
       if match_bond (b,b2) then find_path graph i0 p.id found (i :: rev) bonds else found)

let find_pattern_gen r ids symmetric bond_pattern atom_pattern rep_pattern =
  let found = IntSet.fold ids [] (fun found i ->
    find_path r.graph (if symmetric then i else max_int) i found [] bond_pattern) in
  let found =
    let map = Xlist.fold found StringMap.empty (fun map l -> StringMap.add map (String.concat "-" (Xlist.map l string_of_int)) l) in
    StringMap.fold map [] (fun found _ l -> l :: found) in
  let found =
    if atom_pattern = [] then found else
    Xlist.fold found [] (fun found l ->
      if Xlist.rev_map l (fun i -> (fst r.graph.(i)).name) = atom_pattern then (List.rev l) :: found else found) in
  let found = Xlist.fold found [] (fun found l ->
    try
      let map = Xlist.fold2 l rep_pattern IntMap.empty (fun map i p ->
        IntMap.add_inc map p i (fun j -> if i = j then i else raise Not_found)) in
      let _ = IntMap.fold map IntSet.empty (fun set _ i ->
        if IntSet.mem set i then raise Not_found else IntSet.add set i) in
      l :: found
    with Not_found -> found) in
(*   Xlist.iter found (fun l -> print_endline (String.concat "-" (Xlist.map l string_of_int)));    *)
  found

let find_pattern r symmetric bond_pattern atom_pattern rep_pattern =
  find_pattern_gen r r.reactant_ids symmetric bond_pattern atom_pattern rep_pattern

let add_unbreakable l =
   Collection.sort (Collection.of_list (Xlist.map l (fun (i,j) -> min i j, max i j)))

let find_sigmatropic_kernels r =
(*   print_endline "find_sigmatropic_kernels"; *)
  if Xlist.size r.reactants > 1 then [] else
  let found = find_pattern r true [Double;Single;Single;Single;Double] [] [1;2;3;4;5;6] in
(*  let found = IntSet.fold r.reactant_ids [] (fun found i ->
    find_sigmatropic_kernel r.graph i i found [] [Double;Single;Single;Single;Double]) in
  let found =
    let map = Xlist.fold found StringMap.empty (fun map l -> StringMap.add map (String.concat "-" (Xlist.map l string_of_int)) l) in
    StringMap.fold map [] (fun found _ l -> l :: found) in
  Xlist.iter found (fun l -> print_endline (String.concat "-" (Xlist.map l string_of_int))); *)
  let found = Xlist.fold found [] (fun found l ->
    let n = Xlist.fold l 0 (fun n i -> if IntSet.mem r.cycles i then n+1 else n) in
    if n > 2 then found else l :: found) in
  Xlist.map found (function
    [x1;y1;z1;z2;y2;x2] ->
       let g = add_bonds r.graph [x1,x2,2;x1,y1,-2;x2,y2,-2;y1,z1,2;y2,z2,2;z1,z2,-2] in
       let unb = add_unbreakable [x1,x2;x1,y1;x2,y2;y1,z1;y2,z2] in
       {r with graph=g; broken_bonds=1; msg="Sigmatropic" :: r.msg; unbreakable=unb}
  | _ -> failwith "find_sigmatropic_kernels")

let find_cycloaddition_kernels r =
(*   print_endline "find_cycloaddition_kernels"; *)
  let found = find_pattern r true [Double;Single;Double] [] [1;2;3;4] in
  let found2 = find_pattern r false [Double] [] [1;2] in
  let found = Xlist.rev_map (Xlist.multiply_list [found;found2]) (function [a;b] -> a,b | _ -> failwith "find_cycloaddition_kernels") in
  let found = Xlist.fold found [] (fun found (la,lb) ->
(*     if List.hd la <> 6 then found else (* FIXME!!! *)  *)
(*     if List.hd lb <> 28 && List.hd lb <> 29 then found else (* FIXME!!! *)  *)
    if not (IntSet.is_empty (IntSet.intersection (IntSet.of_list la) (IntSet.of_list lb))) then found else
    (*let n = Xlist.fold la 0 (fun n i -> if IntSet.mem r.cycles i then n+1 else n) in
    if n > 2 then found else*) (la,lb) :: (la,List.rev lb) :: found) in
(*   Xlist.iter found (fun (la,lb) -> print_endline ((String.concat "-" (Xlist.map la string_of_int)) ^ " " ^ (String.concat "-" (Xlist.map lb string_of_int))));  *)
  Xlist.map found (function
    [x1;y1;y2;x2],[z1;z2]->
       let g = add_bonds r.graph [x1,y1,-2;x2,y2,-2;y1,y2,2;x1,z1,2;x2,z2,2;z1,z2,-2] in
       let unb = add_unbreakable [x1,y1;y1,y2;x2,y2;x2,z2;z1,z2;x1,z1] in
(*        let unb_msg = String.concat " " (Xlist.map unb (fun (i,j) -> sprintf "%d-%d" i j)) in *)
       {r with graph=g; broken_bonds=1; msg="Cycloaddition" (*:: unb_msg*) :: r.msg;
           unbreakable=unb}
  | _ -> failwith "find_cycloaddition_kernels")

let find_pummerer_kernels r =
(*   print_endline "find_pummerer_kernels"; *)
  let found = find_pattern r false [Double;Single;Single;Double] ["O";"C";"O";"C";"O"] [1;2;3;4;5] in
  let found2 = find_pattern r false [Double;Single] ["O";"S";"C"] [1;2;3] in
  let found = Xlist.rev_map (Xlist.multiply_list [found;found2]) (function [a;b] -> a,b | _ -> failwith "find_pummerer_kernels") in
  let found = Xlist.fold found [] (fun found (la,lb) ->
    if not (IntSet.is_empty (IntSet.intersection (IntSet.of_list la) (IntSet.of_list lb))) then found else
    let n = Xlist.fold la 0 (fun n i -> if IntSet.mem r.cycles i then n+1 else n) in
    if n > 2 then found else (la,lb) :: found) in
(*   Xlist.iter found (fun (la,lb) -> print_endline ((String.concat "-" (Xlist.map la string_of_int)) ^ " " ^ (String.concat "-" (Xlist.map lb string_of_int))));  *)
  Xlist.map found (function
    [o1;c1;o3;c2;o2],[o4;s1;c3]->
       let g = add_bonds r.graph [c1,o3,-2;s1,o4,-4;o3,c3,2;c1,o4,2] in
       let unb = add_unbreakable [o3,c3;c1,o4;s1,c3] in
       {r with graph=g; broken_bonds=1; msg="Pummerer Rearrangement" :: r.msg; unbreakable=unb}
  | _ -> failwith "find_pummerer_kernels")

let find_wmr_kernels r =
(*   print_endline "find_wmr_kernels"; *)
  let found = find_pattern r false [Double;Single;Single;Single;Double;Single;Single] ["C";"C";"C";"C";"C";"C";"C";"O"] [1;2;3;4;5;6;7;8] in
  let found = Xlist.fold found [] (fun found l ->
    let n = Xlist.fold l 0 (fun n i -> if IntSet.mem r.cycles i then n+1 else n) in
    if n > 2 then found else l :: found) in
  Xlist.map found (function
    [x1;y1;z1;z2;y2;x2;_;_] ->
       let g = add_bonds r.graph [x1,x2,2;x2,y2,-2] in
       let unb = add_unbreakable [x1,y1;y1,z1;z1,z2;z2,y2;y2,x2;x2,x1] in
(*        let unb_msg = String.concat " " (Xlist.map unb (fun (i,j) -> sprintf "%d-%d" i j)) in *)
       {r with graph=g; broken_bonds=1; msg="Wagner-Meerwein Rearrangement" (*:: unb_msg*) :: r.msg; unbreakable=unb}
  | _ -> failwith "find_wmr_kernels")

let find_prins_kernels r =
(*   print_endline "find_prins_kernels"; *)
  let found = find_pattern r false [Double;Single;Single;Single] ["C";"C";"C";"C";"O"] [1;2;3;4;5] in
  let found2 = find_pattern r false [Double] ["C";"O"] [1;2] in (* FIXME: uprościłem wzorzec, miało być "[{C;H}:7]-[C:6](=[O:8])-[{C;H}:9]" *)
  let found = Xlist.rev_map (Xlist.multiply_list [found;found2]) (function [a;b] -> a,b | _ -> failwith "find_prins_kernels") in
  let found = Xlist.fold found [] (fun found (la,lb) ->
    if not (IntSet.is_empty (IntSet.intersection (IntSet.of_list la) (IntSet.of_list lb))) then found else
    let n = Xlist.fold la 0 (fun n i -> if IntSet.mem r.cycles i then n+1 else n) in
    if n > 2 then found else (la,lb) :: (la,List.rev lb) :: found) in
(*   Xlist.iter found (fun (la,lb) -> print_endline ((String.concat "-" (Xlist.map la string_of_int)) ^ " " ^ (String.concat "-" (Xlist.map lb string_of_int))));  *)
  Xlist.map found (function
    [c1;c2;c3;c4;o1],[c5;o2]->
       let g = add_bonds r.graph [c1,c2,-2;c5,o2,-4;o1,c5,2;c1,c5,2] in
       let unb = add_unbreakable [c4,o1;o1,c5] in
       {r with graph=g; broken_bonds=1; msg="Prins cyclization" :: r.msg; unbreakable=unb}
  | _ -> failwith "find_prins_kernels")

let find_prins2_kernels r =
(*   print_endline "find_prins2_kernels"; *)
  let found = find_pattern r false
    [Double;Single;Single;Single;Single;Single;Single]
    ["C";"C";"C";"C";"O";"C";"O";"C"]
    [1;2;3;4;5;6;7;3] in
  Xlist.map found (function
    [c1;c2;c3;c4;o1;c5;o2;_]->
       let g = add_bonds r.graph [c1,c2,-2;c5,o2,-2;c1,c5,2] in
       let unb = add_unbreakable [] in
       {r with graph=g; broken_bonds=1; msg="Prins cyclization 2" :: r.msg; unbreakable=unb}
  | _ -> failwith "find_prins2_kernels")

let find_aldol_kernels r =
(*   print_endline "find_aldol_kernels"; *)
  let found = find_pattern r false [Double;Single;Double] ["O";"C";"C";"O"] [1;2;3;4] in
  let found2 = find_pattern r false [Single;Double;Double;Single] ["C";"C";"O";"C";"C"] [1;2;3;2;4] in
  let found = Xlist.rev_map (Xlist.multiply_list [found;found2]) (function [a;b] -> a,b | _ -> failwith "find_aldol_kernels") in
  let found = Xlist.fold found [] (fun found (la,lb) ->
    if not (IntSet.is_empty (IntSet.intersection (IntSet.of_list la) (IntSet.of_list lb))) then found else
    let n = Xlist.fold lb 0 (fun n i -> if IntSet.mem r.cycles i then n+1 else n) in
    if n > 2 then found else (la,lb) :: found) in
(*   Xlist.iter found (fun (la,lb) -> print_endline ((String.concat "-" (Xlist.map la string_of_int)) ^ " " ^ (String.concat "-" (Xlist.map lb string_of_int))));  *)
  Xlist.map found (function
    [o1;c1;c2;o2],[c3;c4;o5;_;c5]->
       let g = add_bonds r.graph [c1,o1,-4;c2,o2,-4;c1,c3,2;c2,c5,2] in
       {r with graph=g; broken_bonds=1; msg="Aldol Condensation" :: r.msg}
  | _ -> failwith "find_aldol_kernels")

let find_hs_kernels r =
(*   print_endline "find_hs_kernels"; *)
  let found = find_pattern r false [Double;Single;Single] ["C";"C";"C";"Si"] [1;2;3;4] in
  let found2 = find_pattern r false [Double] ["C";"O"] [1;2] in
  let found = Xlist.rev_map (Xlist.multiply_list [found;found2]) (function [a;b] -> a,b | _ -> failwith "find_hs_kernels") in
  let found = Xlist.fold found [] (fun found (la,lb) ->
    if not (IntSet.is_empty (IntSet.intersection (IntSet.of_list la) (IntSet.of_list lb))) then found else (la,lb) :: found) in
(*   Xlist.iter found (fun (la,lb) -> print_endline ((String.concat "-" (Xlist.map la string_of_int)) ^ " " ^ (String.concat "-" (Xlist.map lb string_of_int))));  *)
  Xlist.map found (function
    [c1;c2;c3;si],[c4;o1]->
       let g = add_bonds r.graph [c4,o1,-2;c1,c4,2;c1,c2,-2] in
       {r with graph=g; broken_bonds=1; msg="Hosomi-Sakurai Rearrangement" :: r.msg}
  | _ -> failwith "find_hs_kernels")

let find_wittig_kernels r =
(*   print_endline "find_wittig_kernels"; *)
  let found = find_pattern r false [Single;Single;Single;Single;Double] ["Sn";"C";"O";"C";"C";"C"] [1;2;3;4;5;6] in
(*  let found = Xlist.fold found [] (fun found l ->
    let n = Xlist.fold l 0 (fun n i -> if IntSet.mem r.cycles i then n+1 else n) in
    if n > 2 then found else l :: found) in*)
  Xlist.map found (function
    [sn;c1;o1;c2;c3;c4] ->
       let g = add_bonds r.graph [sn,c1,-2;c3,c4,-2;c1,c4,2] in
       {r with graph=g; broken_bonds=1; msg="2-3-Wittig Rearrangement" :: r.msg}
  | _ -> failwith "find_wittig_kernels")

let find_achmetowicz_kernels r =
(*   print_endline "find_achmetowicz_kernels"; *)
  let found = find_pattern r false [Single;Single;Single;Single;Double;Single;Double] ["O";"C";"C";"O";"C";"C";"C";"C"] [1;2;3;4;5;6;7;3] in
  Xlist.map found (function
    [o1;c1;c2;o2;c3;c4;c5;_] ->
       let g = add_bonds r.graph [o2,c3,-2;o1,c3,2] in
       {r with graph=g; broken_bonds=1; msg="Achmetowicz Rearrangement" :: r.msg}
  | _ -> failwith "find_achmetowicz_kernels")

let find_claisen_kernels r =
(*   print_endline "find_claisen_kernels"; *)
  let found = find_pattern r false [Double;Single;Single;Double] ["C";"C";"O";"C";"O"] [1;2;3;4;5] in
  Xlist.map found (function
    [c1;c2;o1;c3;o2] ->
       let g = add_bonds r.graph [c1,c2,-2;c2,o1,2;o1,c3,-2] in
       {r with graph=g; broken_bonds=1; msg="Claisen Rearrangement" :: r.msg}
  | _ -> failwith "find_claisen_kernels")

let find_claisen2_kernels r =
(*   print_endline "find_claisen_kernels"; *)
  let found = find_pattern r false [Double;Single;Single;Single;Aromatic] ["C";"C";"C";"O";"C";"C"] [1;2;3;4;5;6] in
  Xlist.map found (function
    [c1;c2;c3;o1;c4;c5] ->
       let g = add_bonds r.graph [c1,c2,-2;o1,c3,-2;c1,c5,2;c3,c2,2] in
       {r with graph=g; broken_bonds=1; msg="Claisen2 Rearrangement" :: r.msg}
  | _ -> failwith "find_claisen_kernels")

let find_anhydride_kernels r =
  let found = find_pattern r false [Double;Single;Single;Double] ["O";"C";"O";"C";"O"] [1;2;3;4;5] in
  List.flatten (Xlist.map found (function
    [o1;c1;o2;c2;o3] ->
       let g = add_bonds r.graph [c1,o2,-2] in
       let r = {r with graph=g; broken_bonds=1; msg="Anhydride" :: r.msg} in
       (Xlist.map found (function
           [o1;c1a;o2a;c2;o3] ->
              if o2 = o2a && c1 = c1a then r else
                let g = add_bonds r.graph [c1a,o2a,-2] in
                {r with graph=g}
         | _ -> failwith "find_anhydride_kernels"))
  | _ -> failwith "find_anhydride_kernels"))

let find_anhydride_hydrification_kernels r =
  let found = find_pattern r false [Double;Single;Single;Double] ["O";"C";"O";"C";"O"] [1;2;3;4;5] in
  let found2 = find_pattern r false [Double;Single] ["O";"C";"O"] [1;2;3] in
  List.flatten (Xlist.map found (function
    [o1;c1;o2;c2;o3] ->
       let g = add_bonds r.graph [c1,o2,-2] in
       let r = {r with graph=g; broken_bonds=1; msg="Anhydride Hydrification" :: r.msg} in
       (Xlist.fold found2 [] (fun l -> function
           [o4;c4;o5] ->
              if o4 = o1 || o4 = o2 || o4 = o3 || o5 = o1 || o5 = o2 || o5 = o3 || c4 = c1 || c4 = c2 then l else
                let g = add_bonds r.graph [c4,o5,-2;c1,o5,2] in
                {r with graph=g} :: l
         | _ -> failwith "find_anhydride_hydrification_kernels"))
  | _ -> failwith "find_anhydride_hydrification_kernels"))

let find_olefin_metathesis_kernels r =
(*   print_endline "find_cycloaddition_kernels"; *)
  let found = find_pattern r false [Double] ["C";"C"] [1;2] in
  (* let found1,found2 = Xlist.fold found ([],[]) (fun (found1,found2) -> function
      [c1;c2] -> if Xlist.size (snd r.graph.(c1)) = 1 && Xlist.size (snd r.graph.(c2)) = 1 then (c1,c2) :: found1, found2 else found1, (c1,c2) :: found2
    | _ -> failwith "find_olefin_metathesis_kernels") in *)
  let found = Xlist.fold found [](fun found -> function
      [c1;c2] -> if Xlist.size (snd r.graph.(c1)) = 1 then (c1,c2) :: found else found
    | _ -> failwith "find_olefin_metathesis_kernels") in
  let found = Xlist.rev_map (Xlist.multiply_list [found;found(*found1;found2*)]) (function [a;b] -> a,b | _ -> failwith "find_olefin_metathesis_kernels") in
  Xlist.fold found [] (fun found ((c1,c2),(c3,c4)) ->
    if c1=c3 || c2=c4 || c1=c4 || c2=c3 then found else
    let g = add_bonds r.graph [c1,c2,-4;c3,c4,-4;c1,c3,4;c2,c4,4] in
    let unb = add_unbreakable [c1,c3;c2,c4] in
    {r with graph=g; broken_bonds=1; msg="Olefin Metathesis" :: r.msg;
            unbreakable=unb} :: found)

let find_ketal_kernels r =
(*   print_endline "find_aldol_kernels"; *)
  let found = find_pattern r false [Single;Single] ["O";"C";"O"] [1;2;3] in
  let found2 = find_pattern r false [] ["O"] [1] in
  let found2 = Xlist.fold found2 [] (fun found2 -> function
      [o] -> if Xlist.size (snd r.graph.(o)) = 0 then [o] :: found2 else found2
    | _ -> failwith "find_ketal_kernels") in
  let found = Xlist.rev_map (Xlist.multiply_list [found;found2]) (function [a;b] -> a,b | _ -> failwith "find_ketal_kernels") in
  Xlist.map found (function
    [o1;c1;o2],[o]->
       let g = add_bonds r.graph [c1,o,4;c1,o1,-2;c1,o2,-2] in
       {r with graph=g; broken_bonds=1; msg="Ketal" :: r.msg}
  | _ -> failwith "find_ketal_kernels")

let find_sulphur_kernels r =
(*   print_endline "find_hs_kernels"; *)
  let found = find_pattern r false [Single;Double;Double;Double;Double;Single] ["Cl";"S";"O";"S";"O";"S";"C"] [1;2;3;2;4;2;5] in
  let found2 = find_pattern r false [Single] ["C";"O"] [1;2] in
  let found2 = Xlist.fold found2 [] (fun found2 -> function
      [c;o] -> if Xlist.size (snd r.graph.(o)) = 1 then [c;o] :: found2 else found2
    | _ -> failwith "find_sulphur_kernels") in
  let found = Xlist.rev_map (Xlist.multiply_list [found;found2]) (function [a;b] -> a,b | _ -> failwith "find_hs_kernels") in
  let found = Xlist.fold found [] (fun found (la,lb) ->
    if not (IntSet.is_empty (IntSet.intersection (IntSet.of_list la) (IntSet.of_list lb))) then found else (la,lb) :: found) in
  Xlist.map found (function
    [cl;s;o1;_;o2;_;c1],[c2;o3]->
       let g = add_bonds r.graph [cl,s,-2;o3,s,2;c2,o3,-2] in
       {r with graph=g; broken_bonds=0; msg="Sulphur Rearrangement" :: r.msg}
  | _ -> failwith "find_sulphur_kernels")

let find_anionic_kernels r =
(*   print_endline "find_wittig_kernels"; *)
  let found = find_pattern r false [Single;Single;Single;Single;Double] ["C";"C";"O";"C";"C";"C"] [1;2;3;4;5;6] in
  Xlist.map found (function
    [c;c1;o1;c2;c3;c4] ->
       let g = add_bonds r.graph [o1,c2,-2;c3,c4,-2;c,c4,2] in
       {r with graph=g; broken_bonds=1; msg="Anionic Rearrangement" :: r.msg}
  | _ -> failwith "find_anionic_kernelswittig_kernels")

let find_robinson_kernels r =
  let found = find_pattern r false [Single;Double;Double;Single;Double] ["C";"C";"O";"C";"C";"C"] [1;2;3;2;4;5] in
  let found2 = find_pattern r false [Double;Single] ["O";"C";"C"] [1;2;3] in
  let found = Xlist.rev_map (Xlist.multiply_list [found;found2]) (function [a;b] -> a,b | _ -> failwith "find_robinson_kernels") in
  let found = Xlist.fold found [] (fun found (la,lb) ->
    if not (IntSet.is_empty (IntSet.intersection (IntSet.of_list la) (IntSet.of_list lb))) then found else
    (la,lb) :: found) in
  Xlist.map found (function
    [c4;c3;o1;_;c2;c1],[o2;c5;c6]->
       let g = add_bonds r.graph [o2,c5,-4;c1,c2,-2;c4,c5,2;c1,c6,2] in
       let unb = add_unbreakable [c4,c5;c1,c6] in
(*        let unb_msg = String.concat " " (Xlist.map unb (fun (i,j) -> sprintf "%d-%d" i j)) in *)
       {r with graph=g; broken_bonds=2; msg="Robinson Annulation" (*:: unb_msg*) :: r.msg;
           unbreakable=unb}
  | _ -> failwith "find_robinson_kernels")

let find_carbonyl_kernels r =
  (* print_endline "find_carbonyl_kernels"; *)
  let found = find_pattern r false [Double] ["O";"C"] [1;2] in
  let found2 = find_pattern r false [Single] ["N";"C"] [1;2] in
  let found = Xlist.fold found [] (fun found -> function
      [o;c] ->
        let b = Xlist.fold (snd r.graph.(c)) true (fun b -> function
            Double,p -> if p.id = o then b else false
          | Aromatic,_ -> false
          | _,p -> if p.name = "C" then b else false) in
        if b then [o;c] :: found else found
    | _ -> failwith "find_carbonyl_kernels") in
  let found2 = Xlist.fold found2 [] (fun found2 -> function
      [n;c] ->
        let b = Xlist.fold (snd r.graph.(c)) true (fun b -> function
            Double,_ -> false
          | _ -> b) in
        if Xlist.size (snd r.graph.(n)) = 1 && b then [n;c] :: found2 else found2
    | _ -> failwith "find_carbonyl_kernels") in
  let found = Xlist.rev_map (Xlist.multiply_list [found;found2]) (function [a;b] -> a,b | _ -> failwith "find_carbonyl_kernels") in
  let found = Xlist.fold found [] (fun found (la,lb) ->
    if not (IntSet.is_empty (IntSet.intersection (IntSet.of_list la) (IntSet.of_list lb))) then found else
    (la,lb) :: found) in
  Xlist.map found (function
    [o;c1],[n;c2]->
       let g = add_bonds r.graph [o,c1,-4;n,c1,4] in
       (* let unb = add_unbreakable [c4,c5;c1,c6] in *)
(*        let unb_msg = String.concat " " (Xlist.map unb (fun (i,j) -> sprintf "%d-%d" i j)) in *)
       {r with graph=g; broken_bonds=(-3); msg="Carbonyl" (*:: unb_msg*) :: r.msg;
           (*unbreakable=unb*)}
  | _ -> failwith "find_carbonyl_kernels")

let rec set_aromatic_bond_rec g = function
    x :: y :: l ->
       let p,q = g.(x) in
       g.(x) <- p,Xlist.map q (fun (b,a) -> if a.id = y then (Aromatic : bond),a else b,a);
       let p,q = g.(y) in
       g.(y) <- p,Xlist.map q (fun (b,a) -> if a.id = x then (Aromatic : bond),a else b,a);
       set_aromatic_bond_rec g (y :: l)
  | _ -> ()


let set_aromatic_bond g bonds =
  let g = Array.copy g in
  set_aromatic_bond_rec g bonds;
  g

(*let translate_aromaticity2 r =
  let found = find_pattern_gen r (IntSet.union r.reactant_ids r.product_ids) false [Double;Single;Double;Single;Double;Single;Single;Double;Single;Double;Single]
    ["C";"C";"C";"C";"C";"C";"C";"C";"C";"C";"C";"C"] [1;2;3;4;5;6;1;7;8;9;10;6] in
  let found2 = find_pattern_gen r (IntSet.union r.reactant_ids r.product_ids) false [Single;Double;Single;Double;Single;Double;Single;Double;Single;Double;Single]
    ["C";"C";"C";"C";"C";"C";"C";"C";"C";"C";"C";"C"] [1;2;3;4;5;6;1;7;8;9;10;6] in
  fst (Xlist.fold (found @ found2) (r,IntSet.empty) (fun (r,forbidden) l ->
    if Xlist.fold l false (fun b i -> if IntSet.mem forbidden i then true else false) then r,forbidden else
    let forbidden = Xlist.fold l forbidden IntSet.add in
    {r with graph=set_aromatic_bond r.graph l}, forbidden))*)

let translate_aromaticity1 r =
  let found = find_pattern_gen r (IntSet.union r.reactant_ids r.product_ids) false [Double;Single;Double;Single;Double;Single] ["C";"C";"C";"C";"C";"C";"C"] [1;2;3;4;5;6;1] in
  fst (Xlist.fold found (r,IntSet.empty) (fun (r,forbidden) l ->
    if Xlist.fold l false (fun b i -> if IntSet.mem forbidden i then true else false) then r,forbidden else
    let forbidden = Xlist.fold l forbidden IntSet.add in
    let g = set_aromatic_bond r.graph l in
    {r with graph=g; original_graph=g}, forbidden))

let translate_aromaticity2 r =
  let found = find_pattern_gen r (IntSet.union r.reactant_ids r.product_ids) false [Aromatic;Single;Double;Single;Double;Single] ["C";"C";"C";"C";"C";"C";"C"] [1;2;3;4;5;6;1] in
  fst (Xlist.fold found (r,IntSet.empty) (fun (r,forbidden) l ->
    if Xlist.fold l false (fun b i -> if IntSet.mem forbidden i then true else false) then r,forbidden else
    let forbidden = Xlist.fold l forbidden IntSet.add in
    let g = set_aromatic_bond r.graph l in
    {r with graph=g; original_graph=g}, forbidden))

let translate_aromaticity r =
  let r = translate_aromaticity1 r in
  let r = translate_aromaticity2 r in
  let r = translate_aromaticity2 r in
  let r = translate_aromaticity2 r in
  r

let set_unbreakable r =
  let anhydrides = find_pattern r false [Double;Single;Single;Double] ["O";"C";"O";"C";"O"] [1;2;3;4;5] in (* bezwodniki *)
  let anhydrides = Xlist.fold anhydrides StringSet.empty (fun anhydrides -> function
    [o1;c1;o2;c2;o3] -> StringSet.add anhydrides (Printf.sprintf "%d-%d-%d-%d" o1 c1 o2 c2)
  | _ -> failwith "set_unbreakable") in
(*   print_endline "set_unbreakable";  *)
  let found = find_pattern r false [Single;Single] ["O";"C";"C"] [1;2;3] in
(*   print_endline "set_unbreakable2";  *)
  let l = Xlist.fold found [] (fun l -> function
    [o1;c1;c2] ->
       if Xlist.size (fst r.graph.(o1)).hydrogens = 1 && Xlist.size (fst r.graph.(c1)).hydrogens = 2 then (
(*          print_endline (String.concat "-" (Xlist.map [o1;c1;c2] string_of_int));    *)
         (o1,c1) :: l)
       else l
  | _ -> failwith "set_unbreakable") in
  let found = find_pattern r false [Double;Single;Single] ["O";"C";"O";"C"] [1;2;3;4] in (* FIXME: to łapie też bezwodniki *)
  let l = Xlist.fold found l (fun l -> function
    [o1;c1;o2;c2] -> if StringSet.mem anhydrides (Printf.sprintf "%d-%d-%d-%d" o1 c1 o2 c2) then l else (o2,c2) :: l
  | _ -> failwith "set_unbreakable") in
  {r with unbreakable=add_unbreakable l}

let map_atoms simple_flag re =
      try
        let messages = [] in
(*         let re,messages = Smiles.repair_stoi_reaction re messages in *)
(*        let map = Smiles.calculate_stoi_reaction re in
        let l = List.sort compare (StringMap.fold map [] (fun l name (v,w) ->
          if v = w then l else (name,v,w) :: l)) in
        let messages = if l <> [] then
          let l = Xlist.map l (fun (name,v,w) -> Printf.sprintf "%s: %d$>>$%d" name v w) in
          [Printf.sprintf "INVALID STOICHIOMETRY: %s" (Smiles.escape_string (String.concat " " l))] @ messages
          else messages in
        let re,messages = Smiles.repair_stoi_reaction re messages in*)
(*        let time1 = Sys.time () in
        printf "map_atoms 1\n%!";*)
        let r,messages = Smiles.prepare_record_for_matching re messages in
        let r = translate_aromaticity r in
        let rl = [r] in
(*         let rl = CommonSubstructure.disambiguate_aroma rl in *)
        let rl = if simple_flag then rl else Xlist.map rl set_unbreakable in
(*        let time2 = Sys.time () in
        printf "map_atoms 2 time=%.4f |rl|=%d\n%!" (time2 -. time1) (Xlist.size rl);*)
        let messages,solutions,quality,no_broken_bonds = (*try*) AtomMapping.broken_bonds max_int 6 messages rl (*with _ -> messages,[],max_int,6*) in
        (* let messages,solutions,quality,no_broken_bonds = messages,[],max_int,6 in *)
        let messages,solutions = if simple_flag then messages,solutions else (
(*        let time3 = Sys.time () in
        printf "map_atoms 3 time=%.4f\n%!" (time3 -. time2);*)
(*         let rl = [List.hd rl] in (* FIXME *) *)
          let sigmatropic_rl = (*CommonSubstructure.disambiguate_aroma*) (Xlist.fold rl [] (fun sigmatropic_rl r -> find_sigmatropic_kernels r @ sigmatropic_rl)) in
          let cycloaddition_rl = Xlist.fold rl [] (fun cycloaddition_rl r -> find_cycloaddition_kernels r @ cycloaddition_rl) in
          let pummerer_rl = Xlist.fold rl [] (fun pummerer_rl r -> find_pummerer_kernels r @ pummerer_rl) in
          let wmr_rl = (*CommonSubstructure.disambiguate_aroma*) (Xlist.fold rl [] (fun wmr_rl r -> find_wmr_kernels r @ wmr_rl)) in
          let prins_rl = Xlist.fold rl [] (fun prins_rl r -> find_prins_kernels r @ prins_rl) in
          let prins2_rl = (*CommonSubstructure.disambiguate_aroma*) (Xlist.fold rl [] (fun prins2_rl r -> find_prins2_kernels r @ prins2_rl)) in
          let aldol_rl = Xlist.fold rl [] (fun aldol_rl r -> find_aldol_kernels r @ aldol_rl) in
          let hs_rl = Xlist.fold rl [] (fun hs_rl r -> find_hs_kernels r @ hs_rl) in
          let wittig_rl = Xlist.fold rl [] (fun wittig_rl r -> find_wittig_kernels r @ wittig_rl) in
          let achmetowicz_rl = (*CommonSubstructure.disambiguate_aroma*) (Xlist.fold rl [] (fun achmetowicz_rl r -> find_achmetowicz_kernels r @ achmetowicz_rl)) in
          let claisen_rl = Xlist.fold rl [] (fun claisen_rl r -> find_claisen_kernels r @ claisen_rl) in
          let claisen2_rl = Xlist.fold rl [] (fun claisen_rl r -> find_claisen2_kernels r @ claisen_rl) in
          let anhydride_rl = Xlist.fold rl [] (fun anhydride_rl r -> find_anhydride_kernels r @ anhydride_rl) in
          let anhydride_hydr_rl = Xlist.fold rl [] (fun anhydride_hydr_rl r -> find_anhydride_hydrification_kernels r @ anhydride_hydr_rl) in
          let olefin_metathesis_rl = Xlist.fold rl [] (fun olefin_metathesis_rl r -> find_olefin_metathesis_kernels r @ olefin_metathesis_rl) in
          let ketal_rl = Xlist.fold rl [] (fun ketal_rl r -> find_ketal_kernels r @ ketal_rl) in
          let sulphur_rl = Xlist.fold rl [] (fun ketal_rl r -> find_sulphur_kernels r @ ketal_rl) in
          let anionic_rl = Xlist.fold rl [] (fun ketal_rl r -> find_anionic_kernels r @ ketal_rl) in
          let robinson_rl = Xlist.fold rl [] (fun ketal_rl r -> find_robinson_kernels r @ ketal_rl) in
          let carbonyl_rl = Xlist.fold rl [] (fun ketal_rl r -> find_carbonyl_kernels r @ ketal_rl) in
          let rl(*,messages*) = (*if sigmatropic_rl = [] && cycloaddition_rl = [] && pummerer_rl = [] && wmr_rl = [] && prins_rl = [] && prins2_rl = [] && aldol_rl = [] && hs_rl = [] && wittig_rl = [] && achmetowicz_rl = [] && claisen_rl = [] then*) (*CommonSubstructure.disambiguate_aroma rl*)(*,messages*)
          (*else*)(*@*) sigmatropic_rl @ cycloaddition_rl @ pummerer_rl @ wmr_rl @ prins_rl @ prins2_rl @ aldol_rl @ hs_rl @ wittig_rl @
          achmetowicz_rl @ claisen_rl @ claisen2_rl @ anhydride_rl @ anhydride_hydr_rl @ olefin_metathesis_rl @ ketal_rl @ sulphur_rl @ anionic_rl @ robinson_rl @ carbonyl_rl(*,"Sigmatropic" :: messages*) in
(*         Xlist.iter rl (fun r -> printf "a1 unbreakable=%s\n%!" (String.concat " " (Xlist.map r.unbreakable (fun (i,j) -> sprintf "%d-%d" i j)))); *)
        (*let rl = CommonSubstructure.disambiguate_aroma rl in*) (* FIXME: to przesunąć wyżej *)
(*         Xlist.iter rl (fun r -> printf "a2 unbreakable=%s\n%!" (String.concat " " (Xlist.map r.unbreakable (fun (i,j) -> sprintf "%d-%d" i j)))); *)
(*        let time4a = Sys.time () in
        printf "map_atoms 4a time=%.4f |rl|=%d\n%!" (time4a -. time3) (Xlist.size rl);*)
          let rl = CommonSubstructure.disambiguate_rearrangements rl in
(*        let time4b = Sys.time () in
        printf "map_atoms 4b time=%.4f |rl|=%d\n%!" (time4b -. time4a) (Xlist.size rl);*)
          let messages,solutions2,_,_ = AtomMapping.broken_bonds quality no_broken_bonds messages rl in
(*        let time5 = Sys.time () in
        printf "map_atoms 5 time=%.4f\n%!" (time5 -. time4b);*)
(*        let messages,solutions = Xlist.fold rl (messages,[]) (fun (messages,solutions) r ->
          let m,s = AtomMapping.broken_bonds r in
          m @ messages, s @ solutions) in*)
          let solutions = AtomMapping.select_minimal (solutions @ solutions2) in
          let solutions = ReactionClasses.disambiguate solutions in
          messages,solutions) in
        let messages = if solutions = [] then "No matching found." :: messages else messages in
(*        let time6 = Sys.time () in
        printf "map_atoms 6 time=%.4f\n%!" (time6 -. time5);*)
        messages,solutions
      with e -> [Printexc.to_string e], []

type validation = Validated | NotValidated | NotInValidateSet

let validate xml core_name =
          let no_valid = 0 in
          let valid_xml, no_valid = try Xml.parse_file ("validation/" ^ core_name ^ ".xml"), no_valid+1 with Xml.File_not_found _ -> Xml.PCData "", no_valid in
          let valid_xml2, no_valid = try Xml.parse_file ("validation2/" ^ core_name ^ ".xml"), no_valid+1 with Xml.File_not_found _ -> Xml.PCData "", no_valid in
          let valid_xml3, no_valid = try Xml.parse_file ("validation3/" ^ core_name ^ ".xml"), no_valid+1 with Xml.File_not_found _ -> Xml.PCData "", no_valid in
          let valid_xml4, no_valid = try Xml.parse_file ("validation4/" ^ core_name ^ ".xml"), no_valid+1 with Xml.File_not_found _ -> Xml.PCData "", no_valid in
          let valid_xml5, no_valid = try Xml.parse_file ("validation5/" ^ core_name ^ ".xml"), no_valid+1 with Xml.File_not_found _ -> Xml.PCData "", no_valid in
          let valid_xml6, no_valid = try Xml.parse_file ("validation6/" ^ core_name ^ ".xml"), no_valid+1 with Xml.File_not_found _ -> Xml.PCData "", no_valid in
          if no_valid = 0 then NotInValidateSet else
            if xml = valid_xml || xml = valid_xml2 || xml = valid_xml3 || xml = valid_xml4 then Validated
            else NotValidated


let map_atoms2 simple_flag prefix filename reactions =
  File.file_out filename (fun file ->
    Printf.fprintf file "%s" (Smiles.a0poster_header "a2");
    let _,names = Xlist.fold reactions (1,[]) (fun (n,names) re ->
      Types.time := Sys.time ();
      Printf.printf "%s %s\n%!" re.rxn_id re.reaction_smile;
      try
        Printf.fprintf file "%s {%s}\\\\\n%!" (make_id re.path re.rxn_id) (Smiles.escape_string re.reaction_smile);
        let messages,solutions = map_atoms simple_flag re in
        Xlist.iter messages (fun msg -> Printf.fprintf file "{%s}\\\\\n" (Smiles.escape_string msg));
        Xlist.fold solutions (n,names) (fun (n,names) (r,labels,_,center_bonds,reid) ->
          let added_stoi = Xlist.fold (r.reactants @ r.products) [] (fun added_stoi m ->
            if m.stoi then (IntSet.to_list m.hf_ids) @ (IntSet.to_list m.ids) @ added_stoi else added_stoi) in
          let labels = Labels.add_invisible_list labels added_stoi in
          let xml = Smiles.reaction_to_xml r messages labels center_bonds reid in
          let core_name = prefix ^ make_id re.path re.rxn_id in
          let n,name =
            if Xlist.size solutions = 1 then n, core_name
            else n+1, core_name ^ "_" ^ string_of_int n in
          match validate xml core_name with
            NotInValidateSet ->
            File.file_out ("results/images/" ^ name ^ ".xml") (fun xml_file ->
              Printf.fprintf xml_file "%s\n" (Xml.to_string_fmt xml));
(*           Printf.fprintf file "%s\\\\\n" (Smiles.escape_string (ReactionClasses.to_string (ReactionClasses.classify (r,labels,center_bonds,reid)))); *)
            Xlist.iter r.msg (fun msg -> Printf.fprintf file "{%s}\\\\\n" (Smiles.escape_string msg));
            Printf.fprintf file "\\includegraphics[scale=0.3]{images/%s.png}\\\\\n" name;
            n,(name ^ ".xml") :: names
          | Validated ->
              Printf.fprintf file "Correctly validated\\\\\n";
              n, names
          | NotValidated ->
              Printf.fprintf file "NOT VALIDATED\\\\\n";
              File.file_out ("results/images/" ^ name ^ ".xml") (fun xml_file ->
                Printf.fprintf xml_file "%s\n" (Xml.to_string_fmt xml));
              Xlist.iter r.msg (fun msg -> Printf.fprintf file "{%s}\\\\\n" (Smiles.escape_string msg));
              Printf.fprintf file "\\includegraphics[scale=0.3]{images/%s.png}\\\\\n" name;
              n,(name ^ ".xml") :: names)
      with e -> Printf.fprintf file "{%s}\\\\\n" (Smiles.escape_string (Printexc.to_string e));n,names) in
    Sys.chdir "results/images";
    ignore (Sys.command ("phantomjs mapping-dl.js " ^ String.concat " " names));
    Sys.chdir "../..";
    Printf.fprintf file "%s" Smiles.trailer)

let execution n_workers work (*file*) =
  let size = Xlist.size work in
  let results = Array.make (size+1) (empty_record,[],[]) in
  let r = ref (size + n_workers) in
  let size = string_of_int size in
(*   let output = ref [] in *)
  let _,work = Xlist.fold work (1,[]) (fun (id,work) t ->
    id+1, ((*string_of_int*) id, t) :: work) in
  let work = ref (List.rev work) in
  let id = string_of_int (Unix.getpid ()) in
  let io_list = Int.fold 1 n_workers [] (fun io_list _ ->
    print_endline (id ^ " create_worker");
    let in_chan,out_chan = Unix.open_process "./worker" in
    let descr = Unix.descr_of_in_channel in_chan in
    (in_chan,out_chan,descr) :: io_list) in
  let descr_list = Xlist.map io_list (fun (_,_,descr) -> descr) in
  while !r <> 0 do
    print_endline (id ^ " Unix.select");
    let list,_,_ = Unix.select descr_list [] [] (-1.) in
    print_endline (id ^ " selected " ^ (string_of_int (Xlist.size list)));
    Xlist.iter list (fun descr2 ->
      decr r;
      Xlist.iter io_list (fun (in_chan,out_chan,descr) ->
        if descr = descr2 then (
          let idw = match Marshal.from_channel in_chan with
            Ready_to_work idw ->
              print_endline (idw ^ " ready");
              idw
          | Work_done (idw,(id_work,re,messages,solutions)) ->
              print_endline (idw ^ " work done " ^ re.rxn_id ^ " " ^ re.reaction_smile);
              results.(id_work) <- re,messages,solutions;
(*               output := s :: (!output); *)
(*               Exec.print_result file s; *)
(*               sum_result := Exec.add_result !sum_result s; *)
(*               Exec.print_sum_result file !sum_result; *)
              idw in
          match !work with
            (id,params) :: l ->
              Marshal.to_channel out_chan (Work_with (id,params)) [Marshal.No_sharing];
              flush out_chan;
              print_endline (idw ^ " scheduled " ^ string_of_int id ^ " of " ^ size);
              work := l
          | [] ->
              Marshal.to_channel out_chan Kill_yourself [Marshal.No_sharing];
              print_endline (idw ^ " finished"))))
  done;
  print_endline (id ^ " exit");
  List.rev (Int.fold 1 (Array.length results - 1) [] (fun l i -> results.(i) :: l))


let map_atoms2_distr no_processors simple_flag prefix filename reactions =
  let reactions = List.rev (Xlist.rev_map reactions (fun re -> simple_flag,re)) in
  File.file_out filename (fun file ->
    Printf.fprintf file "%s" (Smiles.a0poster_header "a2");
    let parsed_reactions = execution no_processors reactions in
    let _,names = Xlist.fold parsed_reactions (1,[]) (fun (n,names) (re,messages,solutions) ->
(*       Printf.printf "%s %s\n%!" re.rxn_id re.reaction_smile;  *)
      try
        let b = if Xlist.size solutions = 0 then false else
          Xlist.fold solutions true (fun b (r,labels,_,center_bonds,reid) ->
            let xml = Smiles.reaction_to_xml r messages labels center_bonds reid in
            let core_name = prefix ^ make_id re.path re.rxn_id in
            if validate xml core_name = Validated then b else false) in
        if b then (n,names) else (
        Printf.fprintf file "%s {%s}\\\\\n%!" (make_id re.path re.rxn_id) (Smiles.escape_string re.reaction_smile);
        Xlist.iter messages (fun msg -> Printf.fprintf file "{%s}\\\\\n" (Smiles.escape_string msg));
        let n,names = Xlist.fold solutions (n,names) (fun (n,names) (r,labels,_,center_bonds,reid) ->
          let xml = Smiles.reaction_to_xml r messages labels center_bonds reid in
          let core_name = prefix ^ make_id re.path re.rxn_id in
          let n,name =
            if Xlist.size solutions = 1 then n, core_name
            else n+1, core_name ^ "_" ^ string_of_int n in
          match validate xml core_name with
            NotInValidateSet ->
            File.file_out ("results/images/" ^ name ^ ".xml") (fun xml_file ->
              Printf.fprintf xml_file "%s\n" (Xml.to_string_fmt xml));
(*           Printf.fprintf file "%s\\\\\n" (Smiles.escape_string (ReactionClasses.to_string (ReactionClasses.classify (r,labels,center_bonds,reid)))); *)
            Xlist.iter r.msg (fun msg -> Printf.fprintf file "{%s}\\\\\n" (Smiles.escape_string msg));
            Printf.fprintf file "\\includegraphics[scale=0.3]{images/%s.png}\\\\\n" name;
            n,(name ^ ".xml") :: names
          | Validated ->
              Printf.fprintf file "Correctly validated\\\\\n";
              n, names
          | NotValidated ->
              Printf.fprintf file "NOT VALIDATED\\\\\n";
              File.file_out ("results/images/" ^ name ^ ".xml") (fun xml_file ->
                Printf.fprintf xml_file "%s\n" (Xml.to_string_fmt xml));
              Xlist.iter r.msg (fun msg -> Printf.fprintf file "{%s}\\\\\n" (Smiles.escape_string msg));
              Printf.fprintf file "\\includegraphics[scale=0.3]{images/%s.png}\\\\\n" name;
              n,(name ^ ".xml") :: names) in
          n,names)
      with e -> Printf.fprintf file "{%s}\\\\\n" (Smiles.escape_string (Printexc.to_string e));n,names) in
    Sys.chdir "results/images";
    ignore (Sys.command ("phantomjs mapping-dl.js " ^ String.concat " " names));
    Sys.chdir "../..";
    Printf.fprintf file "%s" Smiles.trailer)

let map_atoms2_distr2 no_processors simple_flag no_reaction no_files_per_dir prefix filename reactions =
  let reactions = List.rev (Xlist.rev_map reactions (fun re -> simple_flag,re)) in
  File.file_out filename (fun file ->
    Printf.fprintf file "%s" (Smiles.a0poster_header "a2");
    let parsed_reactions = execution no_processors reactions in
    let _,_,_ = Xlist.fold parsed_reactions (no_reaction,1,[]) (fun (no_reaction,n,names) (re,messages,solutions) ->
(*       Printf.printf "%s %s\n%!" re.rxn_id re.reaction_smile;  *)
      try
        let b = if Xlist.size solutions = 0 then false else
          Xlist.fold solutions true (fun b (r,labels,_,center_bonds,reid) ->
            let xml = Smiles.reaction_to_xml r messages labels center_bonds reid in
            let core_name = prefix ^ make_id re.path re.rxn_id in
            if validate xml core_name = Validated then b else false) in
        if b then (no_reaction+1,n,names) else (
        Printf.fprintf file "%s {%s}\\\\\n%!" (make_id re.path re.rxn_id) (Smiles.escape_string re.reaction_smile);
        Xlist.iter messages (fun msg -> Printf.fprintf file "{%s}\\\\\n" (Smiles.escape_string msg));
        let dir_name = Printf.sprintf "%s%d" prefix (no_reaction/no_files_per_dir) in
(*         Printf.printf "map_atoms2_distr %s %d %s\n%!" re.rxn_id no_reaction dir_name; *)
        ignore (Sys.command ("mkdir -p results/images/" ^ dir_name));
(*         print_endline "map_atoms2_distr"; *)
        let n,names = Xlist.fold solutions (n,names) (fun (n,names) (r,labels,_,center_bonds,reid) ->
          let xml = Smiles.reaction_to_xml r messages labels center_bonds reid in
          let core_name = prefix ^ make_id re.path re.rxn_id in
          let n,name =
            if Xlist.size solutions = 1 then n, core_name
            else n+1, core_name ^ "_" ^ string_of_int n in
          match validate xml core_name with
            NotInValidateSet ->
            File.file_out ("results/images/" ^ dir_name ^ "/" ^ name ^ ".xml") (fun xml_file ->
              Printf.fprintf xml_file "%s\n" (Xml.to_string_fmt xml));
(*           Printf.fprintf file "%s\\\\\n" (Smiles.escape_string (ReactionClasses.to_string (ReactionClasses.classify (r,labels,center_bonds,reid)))); *)
            Xlist.iter r.msg (fun msg -> Printf.fprintf file "{%s}\\\\\n" (Smiles.escape_string msg));
            Printf.fprintf file "\\includegraphics[scale=0.3]{images/%s/%s.png}\\\\\n" dir_name name;
            n,(name ^ ".xml") :: names
          | Validated ->
              Printf.fprintf file "Correctly validated\\\\\n";
              n, names
          | NotValidated ->
              Printf.fprintf file "NOT VALIDATED\\\\\n";
              File.file_out ("results/images/" ^ dir_name ^ "/" ^ name ^ ".xml") (fun xml_file ->
                Printf.fprintf xml_file "%s\n" (Xml.to_string_fmt xml));
              Xlist.iter r.msg (fun msg -> Printf.fprintf file "{%s}\\\\\n" (Smiles.escape_string msg));
              Printf.fprintf file "\\includegraphics[scale=0.3]{images/%s/%s.png}\\\\\n" dir_name name;
              n,(name ^ ".xml") :: names) in
          no_reaction+1,n,names)
      with e -> Printf.fprintf file "{%s}\\\\\n" (Smiles.escape_string (Printexc.to_string e));no_reaction+1,n,names) in
    Printf.fprintf file "%s" Smiles.trailer)

let map_atoms2_distr3 no_processors simple_flag reactions =
  let reactions = List.rev (Xlist.rev_map reactions (fun re -> simple_flag,re)) in
  let _ = execution no_processors reactions in
  ()



let print_labeled_reactions prefix filename reactions =
  File.file_out filename (fun file ->
    Printf.fprintf file "%s" (Smiles.a0poster_header "a0");
    let _,names = Xlist.fold reactions (1,[]) (fun (n,names) re ->
      Printf.fprintf file "%s {%s}\\\\\n" (make_id re.path re.rxn_id) (Smiles.escape_string re.reaction_smile);
      Types.time := Sys.time ();
      try
      let r,messages = Smiles.prepare_record_for_matching re [] in
      let r = translate_aromaticity r in
      let r = set_unbreakable r in
      let result,history,alt_labels_list = CommonSubstructure.multilevel_label_reaction3 r [] [] in
      let alt_labels_list = if alt_labels_list = [] then [r.empty_labels] else alt_labels_list in
      Xlist.iter messages (fun msg -> Printf.fprintf file "{%s}\\\\\n" (Smiles.escape_string msg));
      Xlist.fold alt_labels_list (n,names) (fun (n,names) labels ->
        let l = Xlist.fold r.reactants [] (fun l m -> set_invisible_hydrogens IntSet.empty l m.tree) in
        let l = Xlist.fold r.products l (fun l m -> set_invisible_hydrogens IntSet.empty l m.tree) in
        let labels = Labels.add_invisible_list labels l in
        let name = prefix ^ string_of_int n in
(*        File.file_out ("results/images/" ^ name ^ ".xml") (fun xml_file ->
          Printf.fprintf xml_file "%s\n" (Xml.to_string_fmt (Smiles.reaction_to_xml r [] labels IntMap.empty IntMap.empty)));*)
(*         Printf.printf "%d %d\n" (IntSet.size r.reactant_ids) (IntSet.size r.product_ids); *)
        File.file_out ("results/images/" ^ name ^ "r.gv") (fun gv_file ->
          Smiles.reaction_to_gv gv_file 1 r.ids_reactants r [] labels IntMap.empty IntMap.empty);
        File.file_out ("results/images/" ^ name ^ "p.gv") (fun gv_file ->
          Smiles.reaction_to_gv gv_file r.ids_reactants r.ids_products r [] labels IntMap.empty IntMap.empty);
        Printf.fprintf file "\\includegraphics[scale=0.3]{images/%sr.png}$\\to$\n" name;
        Printf.fprintf file "\\includegraphics[scale=0.3]{images/%sp.png}\\\\\n" name;
        n+1,(name ^ "p"(* ^ ".xml"*)) :: (name ^ "r") ::  names)
      with e -> Printf.fprintf file "{%s}\\\\\n" (Smiles.escape_string (Printexc.to_string e));n,names) in
    Sys.chdir "results/images";
(*     ignore (Sys.command ("phantomjs mapping-dl.js " ^ String.concat " " names)); *)
(*     Xlist.iter names (fun name -> ignore (Sys.command ("neato -Tpng " ^ name ^ ".gv -o " ^ name ^ ".png"))); *)
    Xlist.iter names (fun name -> ignore (Sys.command ("neato -Tpng " ^ name ^ ".gv -o " ^ name ^ ".png")));
    Sys.chdir "../..";
    Printf.fprintf file "%s" Smiles.trailer)

let print_mapped_reactions prefix filename reactions =
  File.file_out filename (fun file ->
    Printf.fprintf file "%s" (Smiles.a0poster_header "a0");
    let _,names = Xlist.fold reactions (1,[]) (fun (n,names) re ->
      Printf.fprintf file "%s {%s}\\\\\n" (make_id re.path re.rxn_id) (Smiles.escape_string re.reaction_smile);
      Types.time := Sys.time ();
      try
      let r,messages = Smiles.prepare_record_for_matching re [] in
      let r = translate_aromaticity r in
      let r = set_unbreakable r in
      let messages,solutions,quality,no_broken_bonds = try AtomMapping.broken_bonds max_int 6 messages [r]
        with e -> (Printexc.to_string e) :: messages,[],0,0 in
      Xlist.iter messages (fun msg -> Printf.fprintf file "{%s}\\\\\n" (Smiles.escape_string msg));
      let added_stoi = Xlist.fold (r.reactants @ r.products) [] (fun added_stoi m ->
        if m.stoi then (IntSet.to_list m.hf_ids) @ (IntSet.to_list m.ids) @ added_stoi else added_stoi) in
      if solutions = [] then
        let l = Xlist.fold r.reactants added_stoi (fun l m -> set_invisible_hydrogens IntSet.empty l m.tree) in
        let l = Xlist.fold r.products l (fun l m -> set_invisible_hydrogens IntSet.empty l m.tree) in
        let labels = Labels.add_invisible_list r.empty_labels l in
        let name = prefix ^ string_of_int n in
        let reactant_reid = IntSet.fold r.reactant_ids IntMap.empty (fun map i -> IntMap.add map i 0) in
        let product_reid = IntSet.fold r.product_ids IntMap.empty (fun map i -> IntMap.add map i 0) in
        File.file_out ("results/images/" ^ name ^ "r.gv") (fun gv_file ->
          Smiles.reaction_to_gv gv_file 1 r.ids_reactants r [] labels IntMap.empty reactant_reid);
        File.file_out ("results/images/" ^ name ^ "p.gv") (fun gv_file ->
          Smiles.reaction_to_gv gv_file r.ids_reactants r.ids_products r [] labels IntMap.empty product_reid);
        Printf.fprintf file "\\includegraphics[scale=0.3]{images/%sr.png}$\\to$\n" name;
        Printf.fprintf file "\\includegraphics[scale=0.3]{images/%sp.png}\\\\\n" name;
        n+1,(name ^ "p"(* ^ ".xml"*)) :: (name ^ "r") ::  names else
      Xlist.fold solutions (n,names) (fun (n,names) (r,labels,_,center_bonds,reid) ->
        let labels = Labels.add_invisible_list labels added_stoi in
        let name = prefix ^ string_of_int n in
(*        File.file_out ("results/images/" ^ name ^ ".xml") (fun xml_file ->
          Printf.fprintf xml_file "%s\n" (Xml.to_string_fmt (Smiles.reaction_to_xml r [] labels IntMap.empty IntMap.empty)));*)
(*         Printf.printf "%d %d\n" (IntSet.size r.reactant_ids) (IntSet.size r.product_ids); *)
        File.file_out ("results/images/" ^ name ^ "r.gv") (fun gv_file ->
          Smiles.reaction_to_gv gv_file 1 r.ids_reactants r [] labels center_bonds reid);
        File.file_out ("results/images/" ^ name ^ "p.gv") (fun gv_file ->
          Smiles.reaction_to_gv gv_file r.ids_reactants r.ids_products r [] labels center_bonds reid);
        Printf.fprintf file "\\includegraphics[scale=0.3]{images/%sr.png}$\\to$\n" name;
        Printf.fprintf file "\\includegraphics[scale=0.3]{images/%sp.png}\\\\\n" name;
        n+1,(name ^ "p"(* ^ ".xml"*)) :: (name ^ "r") ::  names)
      with e -> Printf.fprintf file "{%s}\\\\\n" (Smiles.escape_string (Printexc.to_string e));n,names) in
    Sys.chdir "results/images";
(*     ignore (Sys.command ("phantomjs mapping-dl.js " ^ String.concat " " names)); *)
(*     Xlist.iter names (fun name -> ignore (Sys.command ("neato -Tpng " ^ name ^ ".gv -o " ^ name ^ ".png"))); *)
    Xlist.iter names (fun name -> ignore (Sys.command ("neato -Tpng " ^ name ^ ".gv -o " ^ name ^ ".png")));
    Sys.chdir "../..";
    Printf.fprintf file "%s" Smiles.trailer)

let print_mapped_reactions_rearr simple_flag prefix filename reactions =
  File.file_out filename (fun file ->
    Printf.fprintf file "%s" (Smiles.a0poster_header "a0");
    let _,names = Xlist.fold reactions (1,[]) (fun (n,names) re ->
      Types.time := Sys.time ();
      let r0,_ = Smiles.prepare_record_for_matching re [] in
      let r0 = {r0 with graph=AtomMapping.expand_hydrogens_and_fluors r0.graph;
                      original_graph=AtomMapping.expand_hydrogens_and_fluors r0.original_graph;
                      reactants=Xlist.map r0.reactants (fun m -> {m with ids = IntSet.union m.ids m.hf_ids});
                      products=Xlist.map r0.products (fun m -> {m with ids = IntSet.union m.ids m.hf_ids});
                      reactant_ids=Xlist.fold r0.reactants IntSet.empty (fun set m -> IntSet.union set m.ids);
                      product_ids=Xlist.fold r0.products IntSet.empty (fun set m -> IntSet.union set m.ids)} in
      Printf.printf "%s %s\n%!" re.rxn_id re.reaction_smile;
      try
        Printf.fprintf file "%s {%s}\\\\\n%!" (make_id re.path re.rxn_id) (Smiles.escape_string re.reaction_smile);
        let messages,solutions = map_atoms simple_flag re in
        Xlist.iter messages (fun msg -> Printf.fprintf file "{%s}\\\\\n" (Smiles.escape_string msg));
        if solutions = [] then
          let r,messages = Smiles.prepare_record_for_matching re [] in
          let r = translate_aromaticity r in
          let added_stoi = Xlist.fold (r.reactants @ r.products) [] (fun added_stoi m ->
            if m.stoi then (IntSet.to_list m.hf_ids) @ (IntSet.to_list m.ids) @ added_stoi else added_stoi) in
          let l = Xlist.fold r.reactants added_stoi (fun l m -> set_invisible_hydrogens IntSet.empty l m.tree) in
          let l = Xlist.fold r.products l (fun l m -> set_invisible_hydrogens IntSet.empty l m.tree) in
          let labels = Labels.add_invisible_list r.empty_labels l in
          let name = prefix ^ string_of_int n in
          let reactant_reid = IntSet.fold r.reactant_ids IntMap.empty (fun map i -> IntMap.add map i 0) in
          let product_reid = IntSet.fold r.product_ids IntMap.empty (fun map i -> IntMap.add map i 0) in
          File.file_out ("results/images/" ^ name ^ "r.gv") (fun gv_file ->
            Smiles.reaction_to_gv gv_file 1 r.ids_reactants r [] labels IntMap.empty reactant_reid);
          File.file_out ("results/images/" ^ name ^ "p.gv") (fun gv_file ->
            Smiles.reaction_to_gv gv_file r.ids_reactants r.ids_products r [] labels IntMap.empty product_reid);
          Printf.fprintf file "\\includegraphics[scale=0.3]{images/%sr.png}$\\to$\n" name;
          Printf.fprintf file "\\includegraphics[scale=0.3]{images/%sp.png}\\\\\n" name;
          n+1,(name ^ "p"(* ^ ".xml"*)) :: (name ^ "r") ::  names else
        Xlist.fold solutions (n,names) (fun (n,names) (r,labels,_,center_bonds,reid) ->
          let added_stoi = Xlist.fold (r.reactants @ r.products) [] (fun added_stoi m ->
            if m.stoi then (IntSet.to_list m.hf_ids) @ (IntSet.to_list m.ids) @ added_stoi else added_stoi) in
          let labels = Labels.add_invisible_list labels added_stoi in
          let name = prefix ^ string_of_int n in
(*        File.file_out ("results/images/" ^ name ^ ".xml") (fun xml_file ->
          Printf.fprintf xml_file "%s\n" (Xml.to_string_fmt (Smiles.reaction_to_xml r [] labels IntMap.empty IntMap.empty)));*)
(*         Printf.printf "%d %d\n" (IntSet.size r.reactant_ids) (IntSet.size r.product_ids); *)
          File.file_out ("results/images/" ^ name ^ "r.gv") (fun gv_file ->
            Smiles.reaction_to_gv gv_file 1 r.ids_reactants r0 [] labels center_bonds reid);
          File.file_out ("results/images/" ^ name ^ "p.gv") (fun gv_file ->
          Smiles.reaction_to_gv gv_file r.ids_reactants r.ids_products r [] labels center_bonds reid);
          Xlist.iter r.msg (fun msg -> Printf.fprintf file "{%s}\\\\\n" (Smiles.escape_string msg));
          Printf.fprintf file "\\includegraphics[scale=0.3]{images/%sr.png}$\\to$\n" name;
          Printf.fprintf file "\\includegraphics[scale=0.3]{images/%sp.png}\\\\\n" name;
          n+1,(name ^ "p"(* ^ ".xml"*)) :: (name ^ "r") ::  names)
      with e -> Printf.fprintf file "{%s}\\\\\n" (Smiles.escape_string (Printexc.to_string e));n,names) in
    Sys.chdir "results/images";
(*     ignore (Sys.command ("phantomjs mapping-dl.js " ^ String.concat " " names)); *)
(*     Xlist.iter names (fun name -> ignore (Sys.command ("neato -Tpng " ^ name ^ ".gv -o " ^ name ^ ".png"))); *)
    Xlist.iter names (fun name -> ignore (Sys.command ("neato -Tpng " ^ name ^ ".gv -o " ^ name ^ ".png")));
    Sys.chdir "../..";
    Printf.fprintf file "%s" Smiles.trailer)

(* cząsteczki we wzorcach nie mogą zaczynać się od R;
   zbiór id atomów ze wzorca musi zawierać wszystkie liczby naturalne od 1 do n *)
let pattern_smarts = [
  "(Trans)Esterification", "[O:2]([R:1])[C:3](=[O:4])[R:5].[O:7]([R:6])[H:8]>>[O:7]([R:6])[C:3](=[O:4])[R:5].[O:2]([R:1])[H:8]",[1;2],[1];
  "(Trans)Esterification", "[O:2]([R:1])[C:3](=[O:4])[R:5].[O:7]([R:6])[H:8]>>[O:7]([R:6])[C:3](=[O:4])[R:5].[O:2]([R:1])[H:8]",[1],[1;2];
  "(Trans)Esterification 2", "[O:2]([R:1])[C:3](=[O:4])[R:5].[O:7]([H:6])[C:8](=[O:9])[R:10]>>[O:7]([H:6])[C:3](=[O:4])[R:5].[O:2]([R:1])[C:8](=[O:9])[R:10]",[1;2],[2];
  "Aminocostam", "[O:1]=[C:2]([O:11][R:12])[C:3]([R:13])=[N:4][N:5]([H:6])[C:7]([R:8])=[N:9][R:10]>>[O:1]=[C:2]1[C:3]([R:13])=[N:4][N:5]=[C:7]([R:8])[N:9]1[R:10].[H:6][O:11][R:12]",[1],[1];
  (* "Sigmatropic Rearrangement", "[C:1]([R:7])([R:8])=[C:2]([R:9])[C:3]([R:10])([R:11])[C:4]([R:12])([R:13])[C:5]([R:14])=[C:6]([R:15])[R:16]>>[C:4]([R:12])([R:13])=[C:5]([R:14])[C:6]([R:15])([R:16])[C:1]([R:7])([R:8])[C:2]([R:9])=[C:3]([R:10])[R:11]",[],[]; występuje postprzegrupowanie *)
  "Amidation", "[O:1]=[C:2]([R:3])[O:4][R:5].[H:6][N:7]([R:8])[R:9]>>[O:1]=[C:2]([R:3])[N:7]([R:8])[R:9].[H:6][O:4][R:5]",[1],[1];
  "Synthesis CBr+HN>>CN+HBr", "[Br:1][C:2]([R:3])([R:4])[R:5].[H:6][N:7]([R:8])[R:9]>>[C:2]([R:3])([R:4])([R:5])[N:7]([R:8])[R:9].[H:6][Br:1]",[1],[1];
  (* "Synthesis HCl", "[Cl:1][R:2].[H:3][R:4]>>[C:2]([R:3])([R:4])([R:5])[N:7]([R:8])[R:9].[H:6][Cl:1]",[],[]; *)
  "Synthesis CCl+HN>>CN+HCl", "[Cl:1][C:2]([R:3])([R:4])[R:5].[H:6][N:7]([R:8])[R:9]>>[C:2]([R:3])([R:4])([R:5])[N:7]([R:8])[R:9].[H:6][Cl:1]",[1],[1];
  (* "Synthesis CCl+HC_2>>CC_2+HCl", "[Cl:1][C:2]([R:3])([R:4])[R:5].[H:6][C:7]([R:8])[R:9]>>[C:2]([R:3])([R:4])([R:5])[C:7]([R:8])[R:9].[H:6][Cl:1]",[],[]; (*FIXME: to nie działa *) *)
  "Synthesis C_2Cl+HO>>CO+HCl", "[Cl:1][C:2]([R:3])[R:5].[H:6][O:7][R:4]>>[C:2]([R:3])([R:5])[O:7][R:4].[H:6][Cl:1]",[1],[1];
  "Synthesis C_2Cl+HN>>CN+HCl", "[Cl:1][C:2]([R:3])[R:5].[H:6][N:7][R:4]>>[C:2]([R:3])([R:5])[N:7][R:4].[H:6][Cl:1]",[1],[1];(*7500*)
  "Synthesis SiCl+HO>>SiO+HCl", "[Cl:1][Si:2]([R:3])([R:4])[R:5].[H:6][O:7][R:8]>>[Si:2]([R:3])([R:4])([R:5])[O:7][R:8].[H:6][Cl:1]",[1],[1];
  "Synthesis CI+HO>>CO+HI", "[I:1][C:2]([R:3])([R:4])[R:5].[H:6][O:7][R:8]>>[C:2]([R:3])([R:4])([R:5])[O:7][R:8].[H:6][I:1]",[1],[1];(*6141*)
  "Ritter 1", "[N:1]#[C:2][R:3].[H:4][O:5][R:6]>>[C:2](=[O:5])([N:1]([H:4])[R:6])[R:3]",[1],[1];
  "Ritter 2", "[N:1]#[C:2][R:3].[O:5]([R:6])[C:7](=[O:8])[R:9].[H:4][O:10][H:11]>>[C:2](=[O:10])([N:1]([H:4])[R:6])[R:3].[O:5]([H:11])[C:7](=[O:8])[R:9]",[1],[1]; (* wersja z etapem wstępnym *)
  (* "Ritter 2", "[N:1]#[C:2][R:3].[O:5]([C:6]([C:12]([H:15])([H:16])[H:17])([C:13]([H:18])([H:19])[H:20])[C:14]([H:21])([H:22])[H:23])[C:7](=[O:8])[R:9].[H:4][O:10][H:11]>>[C:2](=[O:10])([N:1]([H:4])[C:6]([C:12]([H:15])([H:16])[H:17])([C:13]([H:18])([H:19])[H:20])[C:14]([H:21])([H:22])[H:23])[R:3].[O:5]([H:11])[C:7](=[O:8])[R:9]",[1],[1]; (* FIXME: generuje eksplozję kombinatoryczną *)(* wersja z etapem wstępnym *) *)
  (* "Ester Reduction","[O:2]([R:1])[C:3](=[O:4])[R:5].[H:6][H:7].[H:8][H:9]>>[H:9][O:2][R:1].[H:6][O:4][C:3]([H:7])([H:8])[R:5]",[1],[2]; *)
  "Ester Reduction","[O:2]([R:1])[C:3](=[O:4])[R:5].[H+:6].[H-:7].[H-:8].[H+:9]>>[H:9][O:2][R:1].[H:6][O:4][C:3]([H:7])([H:8])[R:5]",[1],[2];
  "Claisen Condensation","[C:1]([R:11])([H:12])([H:10])[C:2](=[O:3])[R:4].[O:5]([R:6])[C:7](=[O:8])[R:9]>>[O:3]=[C:2]([R:4])[C:1]([R:11])([H:12])[C:7](=[O:8])[R:9].[H:10][O:5][R:6]",[1;2],[1];
  "O-Alkilation of Enolates","[O:7]=[C:2]([R:1])[C:3]([H:6])([R:14])[R:4].[S:9](=[O:10])(=[O:11])([O:8][R:5])[O:12][R:13]>>[O:7]([R:13])[C:2]([R:1])=[C:3]([R:14])[R:4].[S:9](=[O:10])(=[O:11])([O:8][R:5])[O:12][H:6]",[2],[1];
  "Amidation 2","[O:1]=[C:2]([R:3])[N:4]([H:5])[H:6].[H:7][N:8]([H:9])[R:10]>>[O:1]=[C:2]([R:3])[N:8]([H:9])[R:10].[H:7][N:4]([H:5])[H:6]",[1],[1];
  "Thioether Synthesis","[O:7]=[C:5]([R:6])[O:4][C:3]([H:1])([H:11])[R:2].[H:10][S:8][R:9]>>[S:8]([R:9])[C:3]([H:1])([H:11])[R:2].[O:7]=[C:5]([R:6])[O:4][H:10]",[1],[1];
  "Condensation with Phosphorus","[P:2](=[O:1])([O:5][R:6])([O:3][R:4])[C:7]([R:8])([R:9])[H:10].[O:11]=[C:12]([R:13])[O:14][R:15]>>[P:2](=[O:1])([O:5][R:6])([O:3][R:4])[C:7]([R:8])([R:9])[C:12](=[O:11])[R:13].[H:10][O:14][R:15]",[1],[1];
  (* "Hydrolysis","[O:2]([R:1])[R:3].[H:4][O:5][H:6]>>[H:4][O:2][R:1].[H:6][O:5][R:3]",[1],[]; *)
  ]

let patterns = Xlist.map pattern_smarts Smiles.prepare_pattern

let print_molecule m =
  Printf.printf "%s%s %s {%s} {%s}\n" (if m.obligatory then "* " else "") m.smiles
    (Smiles.string_of_tree_std_id (Smiles.smile_tree_of_atom_tree m.tree))
    (String.concat "," (Xlist.rev_map (IntSet.to_list m.ids) string_of_int))
    (String.concat "," (Xlist.rev_map (IntSet.to_list m.hf_ids) string_of_int))

let string_of_charge = function
    0 -> ""
  | 1 -> "+"
  | -1 -> "-"
  | n -> if n > 0 then "+" ^ string_of_int n else "-" ^ string_of_int n

let string_of_atom_props p =
  Printf.sprintf "%2d %s%s {%s} {%s}" p.id p.name (string_of_charge p.charge)
    (String.concat "," (Xlist.map p.hydrogens string_of_int))
    (String.concat "," (Xlist.map p.fluors string_of_int))

let print_graph g =
  Array.iteri (fun i (p,q) ->
    Printf.printf "%s %s\n" (string_of_atom_props p)
      (String.concat " " (Xlist.map q (fun (b,p) -> Smiles.string_of_bond b ^ string_of_int p.id)))) g

let print_pattern p =
  Printf.printf "%s\nReactants:\n" p.pat_name;
  Xlist.iter p.reactants print_molecule;
  (* print_graph p.reactant_graph; *)
  Printf.printf "Products:\n";
  Xlist.iter p.products print_molecule;
  (* print_graph p.product_graph; *)
  ()

let print_reaction r =
  Printf.printf "%s\nReactants:\n" r.record.reaction_smile;
  Xlist.iter r.reactants print_molecule;
  (* print_graph p.reactant_graph; *)
  Printf.printf "Products:\n";
  Xlist.iter r.products print_molecule;
  (* print_graph p.product_graph; *)
  Xlist.iter r.msg print_endline;
  ()

let print_matching_arrays reactant_array product_array =
  let n = Array.length reactant_array - 1 in
  Printf.printf " ";
  Int.iter 1 n (fun i -> Printf.printf " %3d" i);
  Printf.printf "\nr";
  Int.iter 1 n (fun i -> Printf.printf " %3d" reactant_array.(i));
  Printf.printf "\np";
  Int.iter 1 n (fun i -> Printf.printf " %3d" product_array.(i));
  Printf.printf "\n"

let rec find_molecule x = function
    [] -> failwith "find_molecule"
  | m :: l -> if IntSet.mem m.ids x then m else find_molecule x l

let create_link_id s x y =
  s ^ string_of_int (min x y) ^ "-" ^ string_of_int (max x y)

let rec create_tree_of_graph graph visited super_id id =
  (*Printf.printf "create_tree_of_graph: super_id=%d id=%d\n" super_id id;
  print_graph graph;
  print_endline (String.concat " " (Xlist.map (IntSet.to_string visited) string_of_int));*)
  if IntSet.mem visited id then Link(create_link_id "G" super_id id),visited else
  let p,l = graph.(id) in
  let visited = IntSet.add visited id in
  let l,visited = Xlist.fold l ([],visited) (fun (l,visited) (b,q) ->
    if q.id = super_id then l,visited else
    let a,visited = create_tree_of_graph graph visited id q.id in
    (b,a) :: l,visited) in
  Atom(p, List.rev l),visited

let rec create_tree graph visited super_id matching = function
    Atom(p,l) ->
      let id = matching.(p.id) in
      if p.name = "R" then
        if l <> [] then failwith "create_tree" else
        if id = 0 then Atom({p with id=(-p.id)},[]),visited
        else create_tree_of_graph graph visited super_id id
      else
        let p = if id = 0 then {p with id=(-p.id)} else fst graph.(id) in
        let l,visited = Xlist.fold l ([],visited) (fun (l,visited) (b,q) ->
          let a,visited = create_tree graph visited id matching q in
          (b,a) :: l,visited) in
        Atom(p, List.rev l),visited
  | Link n -> Link n,visited

let rec manage_links next_link link_map = function
    Atom(p,l) ->
      let l,next_link,link_map = Xlist.fold l ([],next_link,link_map) (fun (l,next_link,link_map) (b,q) ->
          let a,next_link,link_map = manage_links next_link link_map q in
          (b,a) :: l,next_link,link_map) in
      Atom(p, List.rev l),next_link,link_map
  | Link s ->
      if StringMap.mem link_map s then Link(StringMap.find link_map s),next_link,link_map else
      let n = string_of_int next_link in
      Link n, next_link+1, StringMap.add link_map s n

let is_improper = function
    Atom(p,l) ->
      let links = Xlist.fold l StringMap.empty (fun links a -> Smiles.make_links_map links p a) in
      StringMap.fold links false (fun b _ -> function
              (_,_,[_]) -> b
            | _ -> true)
  | Link _ -> true

exception Improper_molecule

let completize_pattern graph p_graph molecules pat straight_matching opposite_matching =
  (* print_endline ("completize_pattern " ^ pat.smiles); *)
  let x = IntSet.min_elt pat.ids in
  let b = straight_matching.(x) in
  (* Printf.printf "x=%d b=%d\n" x b; *)
  if b > 0 then find_molecule b molecules, [] else (
  let visited = Int.fold 1 (Array.length opposite_matching - 1) IntSet.empty (fun visited id ->
    if (fst p_graph.(id)).name = "R" then visited else IntSet.add visited opposite_matching.(id)) in
  let tree,_ = create_tree graph visited 0 opposite_matching pat.tree in
  (* print_endline "completize_pattern 1";
  print_endline (Smiles.string_of_tree_std_id (Smiles.smile_tree_of_atom_tree tree)); *)
  let tree,_,_ = manage_links 1 StringMap.empty tree in
  if is_improper tree then raise Improper_molecule else (* to ma miejsce, gdy reszty mające trafić do różnych molekuł są ze sobą połączone *)
  let smiles = Smiles.string_of_tree_no_id (Smiles.smile_tree_of_atom_tree tree) in
  let mol = {smiles=smiles; smile_tree=SLink ""; tree=tree;
             ids=IntSet.empty; hf_ids=IntSet.empty; stoi=false; obligatory=false; cycles=[]} in
  (* print_endline "synthetized component";
  print_molecule mol; *)
  mol, [smiles])

let rec reid_tree next_id reid_map = function
    Atom(p,l) ->
      if IntMap.mem reid_map p.id then failwith "reid_tree" else
      let reid_map = IntMap.add reid_map p.id next_id in
      let p = {p with id=next_id} in
      let next_id,reid_map,l = Xlist.fold l (next_id+1,reid_map,[]) (fun (next_id,reid_map,l) (b,q) ->
        let next_id,reid_map,a = reid_tree next_id reid_map q in
        next_id,reid_map,(b,a) :: l) in
      next_id,reid_map,Atom(p, List.rev l)
  | Link n -> next_id,reid_map,Link n

let completize r p reactant_matching product_matching =
  (* print_endline "completize";
  print_matching_arrays reactant_matching product_matching; *)
  let reactants, added_reactant_smiles = Xlist.fold p.reactants ([],[]) (fun (reactants,added_reactant_smiles) pat ->
    let mol,smile = completize_pattern r.graph p.reactant_graph r.reactants pat reactant_matching product_matching in
    mol :: reactants, smile @ added_reactant_smiles) in
  let products, added_product_smiles = Xlist.fold p.products ([],[]) (fun (products,added_product_smiles) pat ->
    let mol,smile = completize_pattern r.graph p.product_graph r.products pat product_matching reactant_matching in
    mol :: products, smile @ added_product_smiles) in
  let ids_reactants,reid_reactants,reactants = Xlist.fold reactants (1,IntMap.empty,[]) (fun (next_id,reid_map,reactants) molecule ->
    let next_id,reid_map,tree = reid_tree next_id reid_map molecule.tree in
    next_id,reid_map,{molecule with tree=tree; ids=Smiles.get_tree_ids IntSet.empty tree} :: reactants) in
  let ids_products,reid_products,products = Xlist.fold products (ids_reactants,IntMap.empty,[]) (fun (next_id,reid_map,products) molecule ->
    let next_id,reid_map,tree = reid_tree next_id reid_map molecule.tree in
    next_id,reid_map,{molecule with tree=tree; ids=Smiles.get_tree_ids IntSet.empty tree} :: products) in
  let reactant_reid_matching = Array.make p.ids_reactants 0 in
  let product_reid_matching = Array.make p.ids_reactants 0 in
  let graph = Smiles.make_atom_graph ids_products (Xlist.map (reactants @ products) (fun m -> m.tree)) in
  let reactant_ids = Xlist.fold reactants IntSet.empty (fun set m -> IntSet.union set m.ids) in
  let product_ids = Xlist.fold products IntSet.empty (fun set m -> IntSet.union set m.ids) in
  let empty_labels = Labels.make ids_products in
  let reactant_msg = if added_reactant_smiles = [] then [] else ["Added reactants: " ^ String.concat " " added_reactant_smiles] in
  let product_msg = if added_product_smiles = [] then [] else ["Added products: " ^ String.concat " " added_product_smiles] in
  let reaction = {empty_reaction with
     reactants = reactants;
     products = products;
     ids_reactants = ids_reactants;
     ids_products = ids_products;
     reactant_ids = reactant_ids;
     product_ids = product_ids;
     graph = graph;
     original_graph = graph;
     empty_labels = empty_labels;
     reaction_size = ids_products;
     msg = ["INVALID STOI"] @ reactant_msg @ product_msg @ r.msg;
     record = r.record
     } in
  (* zmiana indeksów uwzględniająca 3 przypadki: zmaczowane, jeden brakujący,  dwa brakujące *)
  (* print_endline "completize 3"; *)
  Int.iter 1 (p.ids_reactants - 1) (fun id ->
    (* print_endline ("completize 2a " ^ string_of_int id); *)
    reactant_reid_matching.(id) <-
      if reactant_matching.(id) = 0 && product_matching.(id) = 0 then IntMap.find reid_reactants (-id) else
      if reactant_matching.(id) = 0 then IntMap.find reid_reactants product_matching.(id) else
      IntMap.find reid_reactants reactant_matching.(id);
    (* print_endline "completize 2b"; *)
    product_reid_matching.(id) <-
      if reactant_matching.(id) = 0 && product_matching.(id) = 0 then IntMap.find reid_products (-id) else
      if product_matching.(id) = 0 then IntMap.find reid_products reactant_matching.(id) else
      IntMap.find reid_products product_matching.(id));
  (* print_endline "reid matching";
  print_matching_arrays reactant_reid_matching product_reid_matching; *)
  reaction,reactant_reid_matching,product_reid_matching

let map_one_side r_hash p_hash r_graph p_graph r_selection p_selection =
  let obl,p_selection = Xlist.fold p_selection ([],[]) (fun (obl,p_selection) m ->
    if m.obligatory then m :: obl,p_selection else obl,m :: p_selection) in
  let obl = Collection.of_list obl in
  let cands = Collection.generate_combinations (Xlist.size r_selection - Collection.size obl) (Collection.of_list p_selection) in
  let cands = Collection.map cands (fun pats -> Collection.sum obl pats) in
  (* Collection.print_cc cands "cands1" "" (fun m -> m.smiles); *)
  let cands = Collection.flatten_map cands (fun pats ->
    Collection.generate_permutations (pats, Collection.of_list r_selection)) in
  (* Collection.print_cc cands "cands2" "cand" (fun (mol,pat) -> mol.smiles ^ " <---> " ^ pat.smiles); *)
  let matching = Collection.flatten_map cands (fun pat_mols ->
    let matching = Collection.map pat_mols (fun (pat,mol) ->
      let c = (*try*) Isomorphism.find_isomorphisms true 0 p_hash r_hash p_graph r_graph
        (Collection.of_list (IntSet.to_list pat.ids))
        (Collection.of_list (IntSet.to_list mol.ids)) (*with e -> (print_endline (Printexc.to_string e); Collection.empty)*) in
        (* Collection.print_cc c "matching1" "" (fun (i,j) -> string_of_int i ^ "-" ^ string_of_int j); *)
      c) in
      (* Collection.print_ccc matching "matching2" "pat-mol" "" (fun (i,j) -> string_of_int i ^ "-" ^ string_of_int j); *)
    if Collection.int_log (Collection.number_of_products matching) > 6 then failwith "map_patterns"(*raise (Ambiguity "too many matchings in map_patterns")*); (* FIXME *)
    let matching = Collection.flatten_product matching in
    (* Collection.print_cc matching "matching3" "" (fun (i,j) -> string_of_int i ^ "-" ^ string_of_int j); *)
    matching) in
  matching

let is_incomplete a =
  Int.fold 1 (Array.length a - 1) false (fun b id -> a.(id) = 0 || b)

let rec count_obligatory = function
    [] -> 0
  | m :: l -> (if m.obligatory then 1 else 0) + count_obligatory l

let map_patterns messages patterns r =
  (* print_reaction r; *)
  let hash = Hash.make_hash Hash.hash_key 0 r.graph r.empty_labels in
  Xlist.fold patterns (messages,[]) (fun (messages,solutions) p ->
    (* print_pattern p; *)
    if Xlist.size r.reactants > Xlist.size p.reactants ||
       Xlist.size r.products > Xlist.size p.products then messages,solutions else
    if Xlist.size r.reactants < count_obligatory p.reactants ||
       Xlist.size r.products < count_obligatory p.products then messages,solutions else (
    let reactant_hash = Hash.make_hash Hash.hash_key 0 p.reactant_graph p.empty_labels in
    let reactant_matchings = map_one_side hash reactant_hash r.graph p.reactant_graph r.reactants p.reactants in
    (* Collection.print_cc reactant_matchings "reactant_matching" "" (fun (i,j) -> string_of_int i ^ "-" ^ string_of_int j); *)
    let product_hash = Hash.make_hash Hash.hash_key 0 p.product_graph p.empty_labels in
    let product_matchings = map_one_side hash product_hash r.graph p.product_graph r.products p.products in
    (* Collection.print_cc product_matchings "product_matching" "" (fun (i,j) -> string_of_int i ^ "-" ^ string_of_int j); *)
    let alt_matchings = Collection.flatten_map reactant_matchings (fun reactant_matching ->
      let reactant_array = Array.make p.ids_reactants 0 in
      Xlist.iter (Collection.to_list reactant_matching) (fun (p_id,r_id) -> reactant_array.(p_id) <- r_id);
      Collection.flatten_map product_matchings (fun product_matching ->
        let product_array = Array.make p.ids_reactants 0 in
        Xlist.iter (Collection.to_list product_matching) (fun (p_id,r_id) -> product_array.(p_id) <- r_id);
        (* print_endline "\nmatching";
        print_matching_arrays reactant_array product_array; *)
        let b = Int.fold 1 (p.ids_reactants - 1) true (fun b id ->
          let reactant_id = reactant_array.(id) in
          let product_id = product_array.(id) in
          if reactant_id = 0 || product_id = 0 then b else
          ((fst r.graph.(reactant_id)).name = (fst r.graph.(product_id)).name) && b) in
        if not b then Collection.empty else
        try
          let incomplete = is_incomplete reactant_array || is_incomplete product_array in
          let reaction,reactant_matching,product_matching =
            if incomplete then completize r p reactant_array product_array
            else r,reactant_array,product_array in
          let matching = Int.fold 1 (p.ids_reactants-1) [] (fun l i ->
            (reactant_matching.(i), product_matching.(i)) :: l) in
          let reaction = {reaction with msg = p.pat_name :: reaction.msg} in
          Collection.singleton (reaction,matching)
        with Improper_molecule -> Collection.empty)) in
    (* print_endline "second stage"; *)
    let alt_matchings = Collection.map alt_matchings (fun (r,matching) ->
      (* print_reaction r; *)
      let graph = Smiles.collapse_hydrogens_and_fluors r.graph in
      let reactants = Xlist.map r.reactants (fun m -> let ids,hf_ids = Smiles.select_non_hf_atoms graph m.ids in {m with ids=ids; hf_ids=hf_ids}) in
      let products = Xlist.map r.products (fun m -> let ids,hf_ids = Smiles.select_non_hf_atoms graph m.ids in {m with ids=ids; hf_ids=hf_ids}) in
      let reactant_ids = Xlist.fold reactants IntSet.empty (fun set m -> IntSet.union set m.ids) in
      let product_ids = Xlist.fold products IntSet.empty (fun set m -> IntSet.union set m.ids) in
      let r = {r with
        reactants = reactants;
        products = products;
        reactant_ids = reactant_ids;
        product_ids = product_ids;
        graph = graph;
        original_graph = graph;
        } in
      r,matching) in
    let messages,solutions2 = AtomMapping.broken_bonds_pat messages alt_matchings in
    messages,solutions2 @ solutions))

let find_aux_molecules l =
  let pri,aux = Xlist.fold l ([],[]) (fun (pri,aux) s ->
    if s = "Cl" then pri, s :: aux else s :: pri, aux) in
  List.rev pri, List.rev aux

let remove_aux_molecules re =
  let reactant_smiles,product_smiles = Smiles.parse_smile_reaction re.reaction_smile in
  let reactant_smiles,aux_reactants = find_aux_molecules reactant_smiles in
  let product_smiles,aux_products = find_aux_molecules product_smiles in
  let reaction_smiles = String.concat "." reactant_smiles ^ ">>" ^ String.concat "." product_smiles in
  let messages =
    (if aux_reactants = [] then [] else ["Removed reactants: " ^ String.concat " " aux_reactants]) @
    (if aux_products = [] then [] else ["Removed products: " ^ String.concat " " aux_products]) in
  let flag = aux_reactants <> [] || aux_products <> [] in
  messages, {re with reaction_smile=reaction_smiles}, flag

let map_atoms_pat re =
  (* print_endline "map_atoms_pat"; *)
       try
        let messages = [] in
        let r,messages = Smiles.prepare_record_for_matching_no_stoi re messages in
        let r = translate_aromaticity r in
        let messages,solutions = map_patterns messages patterns r in
        let messages,solutions = if solutions <> [] then messages,solutions else
          let messages2,re,flag = remove_aux_molecules re in
          if not flag then messages,solutions else
            let r,messages = Smiles.prepare_record_for_matching_no_stoi re messages2 in
            let r = translate_aromaticity r in
            map_patterns messages patterns r in
        let messages = if solutions = [] then "No matching found." :: messages else messages in
        messages,solutions
       with e -> [Printexc.to_string e], []


let print_mapped_reactions_rearr2 prefix filename reactions =
  File.file_out filename (fun file ->
    Printf.fprintf file "%s" (Smiles.a0poster_header "a0");
    let _,names = Xlist.fold reactions (1,[]) (fun (n,names) re ->
      Types.time := Sys.time ();
      (*let r0,_ = Smiles.prepare_record_for_matching re [] in
      let r0 = {r0 with graph=AtomMapping.expand_hydrogens_and_fluors r0.graph;
                      original_graph=AtomMapping.expand_hydrogens_and_fluors r0.original_graph;
                      reactants=Xlist.map r0.reactants (fun m -> {m with ids = IntSet.union m.ids m.hf_ids});
                      products=Xlist.map r0.products (fun m -> {m with ids = IntSet.union m.ids m.hf_ids});
                      reactant_ids=Xlist.fold r0.reactants IntSet.empty (fun set m -> IntSet.union set m.ids);
                      product_ids=Xlist.fold r0.products IntSet.empty (fun set m -> IntSet.union set m.ids)} in*)
      Printf.printf "%s %s\n%!" re.rxn_id re.reaction_smile;
      try
        Printf.fprintf file "%s {%s}\\\\\n%!" (make_id re.path re.rxn_id) (Smiles.escape_string re.reaction_smile);
        let messages,solutions = map_atoms_pat re in
        Xlist.iter messages (fun msg -> Printf.fprintf file "{%s}\\\\\n" (Smiles.escape_string msg));
        (*if solutions = [] then
          let r,messages = Smiles.prepare_record_for_matching re [] in
          let r = translate_aromaticity r in
          let added_stoi = Xlist.fold (r.reactants @ r.products) [] (fun added_stoi m ->
            if m.stoi then (IntSet.to_list m.hf_ids) @ (IntSet.to_list m.ids) @ added_stoi else added_stoi) in
          let l = Xlist.fold r.reactants added_stoi (fun l m -> set_invisible_hydrogens IntSet.empty l m.tree) in
          let l = Xlist.fold r.products l (fun l m -> set_invisible_hydrogens IntSet.empty l m.tree) in
          let labels = Labels.add_invisible_list r.empty_labels l in
          let name = prefix ^ string_of_int n in
          let reactant_reid = IntSet.fold r.reactant_ids IntMap.empty (fun map i -> IntMap.add map i 0) in
          let product_reid = IntSet.fold r.product_ids IntMap.empty (fun map i -> IntMap.add map i 0) in
          File.file_out ("results/images/" ^ name ^ "r.gv") (fun gv_file ->
            Smiles.reaction_to_gv gv_file 1 r.ids_reactants r [] labels IntMap.empty reactant_reid);
          File.file_out ("results/images/" ^ name ^ "p.gv") (fun gv_file ->
            Smiles.reaction_to_gv gv_file r.ids_reactants r.ids_products r [] labels IntMap.empty product_reid);
          Printf.fprintf file "\\includegraphics[scale=0.3]{images/%sr.png}$\\to$\n" name;
          Printf.fprintf file "\\includegraphics[scale=0.3]{images/%sp.png}\\\\\n" name;
          n+1,(name ^ "p"(* ^ ".xml"*)) :: (name ^ "r") ::  names else*)
        Xlist.fold solutions (n,names) (fun (n,names) (r,labels,_,center_bonds,reid) ->
          let added_stoi = Xlist.fold (r.reactants @ r.products) [] (fun added_stoi m ->
            if m.stoi then (IntSet.to_list m.hf_ids) @ (IntSet.to_list m.ids) @ added_stoi else added_stoi) in
          let labels = Labels.add_invisible_list labels added_stoi in
          let name = prefix ^ string_of_int n in
(*        File.file_out ("results/images/" ^ name ^ ".xml") (fun xml_file ->
          Printf.fprintf xml_file "%s\n" (Xml.to_string_fmt (Smiles.reaction_to_xml r [] labels IntMap.empty IntMap.empty)));*)
(*         Printf.printf "%d %d\n" (IntSet.size r.reactant_ids) (IntSet.size r.product_ids); *)
          File.file_out ("results/images/" ^ name ^ "r.gv") (fun gv_file ->
            Smiles.reaction_to_gv gv_file 1 r.ids_reactants r [] labels center_bonds reid);
          File.file_out ("results/images/" ^ name ^ "p.gv") (fun gv_file ->
            Smiles.reaction_to_gv gv_file r.ids_reactants r.ids_products r [] labels center_bonds reid);
          Xlist.iter r.msg (fun msg -> Printf.fprintf file "{%s}\\\\\n" (Smiles.escape_string msg));
          Printf.fprintf file "\\includegraphics[scale=0.3]{images/%sr.png}$\\to$\n" name;
          Printf.fprintf file "\\includegraphics[scale=0.3]{images/%sp.png}\\\\\n" name;
          n+1,(name ^ "p"(* ^ ".xml"*)) :: (name ^ "r") ::  names)
      with e -> (Printf.fprintf file "{%s}\\\\\n" (Smiles.escape_string (Printexc.to_string e));n,names)) in
    Sys.chdir "results/images";
(*     ignore (Sys.command ("phantomjs mapping-dl.js " ^ String.concat " " names)); *)
(*     Xlist.iter names (fun name -> ignore (Sys.command ("neato -Tpng " ^ name ^ ".gv -o " ^ name ^ ".png"))); *)
    Xlist.iter names (fun name -> ignore (Sys.command ("neato -Tpng " ^ name ^ ".gv -o " ^ name ^ ".png")));
    Sys.chdir "../..";
    Printf.fprintf file "%s" Smiles.trailer)
