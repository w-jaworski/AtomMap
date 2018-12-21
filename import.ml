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

let get_host_name () =
  let chan = Unix.open_process_in "uname -n" in
  input_line chan

let nexus_path = "/home/yacheu/Dokumenty/Badania/Chemia/"
let yoshimune_path = "/home/wjaworski/Chemia/"
let terra_path = "/home/yacheu/Dokumenty/Chemia/"
let students_path = "/home/wjaworski/Chemia/"
let wloczykij_path = "/home/wjaworski/Dokumenty/Chemia/"
let kolos_path = "/home/mapper/mapper/"

let current_path =
  match get_host_name () with
    "nexus" -> nexus_path
  | "yoshimune" -> yoshimune_path
  | "terra" -> terra_path
  | "students" -> students_path
  | "wloczykij" -> wloczykij_path
  | "kolos" -> kolos_path
  | "TERRA" -> nexus_path
  | s -> print_endline ("unknown host: " ^ s); ""

let rxn_filename = current_path ^ "src_data/rxn_db.dat"
let rxn_stoi_filename = current_path ^ "src_data/rxn_db_stoi.dat"
(*let paths_filename = current_path ^ "src_data/Mapping_algorithm_train_set1.dat"
let organic_syntheses_filename = current_path ^ "src_data/OrganicSynthesesReactions.dat"
let organic_syntheses2_filename = current_path ^ "src_data/200_OrganicSyntheses_Apr7_2015.dat"
let przyklady_trudne_filename = current_path ^ "src_data/przyklady_trudne.dat"
let przyklady_trudne2_filename = current_path ^ "src_data/reakcje_do_zmaczowania_trudne_2.dat"*)
let paths_filename = "../src_data/Mapping_algorithm_train_set1.dat"
let organic_syntheses_filename = "../src_data/OrganicSynthesesReactions.dat"
let organic_syntheses2_filename = "../src_data/200_OrganicSyntheses_Apr7_2015.dat"
let przyklady_trudne_filename = "../src_data/przyklady_trudne.dat"
let przyklady_trudne2_filename = "../src_data/reakcje_do_zmaczowania_trudne_2.dat"
let losowe_patenty_filename = "../src_data/losowe_patenty.dat"
let rt_filename = current_path ^ "src_data/RT.txt"
let metacyc_dir = current_path ^ "src_data/00_input/"
let w_5_ohinprod_filename = "../src_data/W_5_OHinProd.txt"
let w_5_ohinsub_filename = "../src_data/W_5_OHinSub.txt"
let w_6_ohinprod_filename = "../src_data/W_6_OHinProd.txt"
let w_6_ohinsub_filename = "../src_data/W_6_OHinSub.txt"
let wszystkie_filename = "../wszystkie.csv"

let rxn_molecules_filename = current_path ^ "inf_data/rxn_db_molecules.dat"
let paths_molecules_filename = current_path ^ "inf_data/Mapping_algorithm_train_set1_molecules.dat"
let organic_syntheses_molecules_filename = current_path ^ "inf_data/OrganicSynthesesReactions_molecules.dat"
let organic_syntheses2_molecules_filename = current_path ^ "inf_data/200_OrganicSyntheses_Apr7_2015_molecules.dat"
let example_molecules_filename = "results/molecules2.dat"

(*let rxn_molecules_visualization_filename = current_path ^ "inf_data/rxn_db_molecules_visualization.dat"
let paths_molecules_visualization_filename = current_path ^ "inf_data/Mapping_algorithm_train_set1_molecules_visualization.dat"
let organic_syntheses_molecules_visualization_filename = current_path ^ "inf_data/OrganicSynthesesReactions_molecules_visualization.dat"
let organic_syntheses2_molecules_visualization_filename = current_path ^ "inf_data/200_OrganicSyntheses_Apr7_2015_molecules_visualization.dat"*)

let molecules_visualization_filename = current_path ^ "inf_data/molecules_visualization.dat"

(* let simple_reaction_ids_filename = current_path ^ "inf_data/simple_reactions.dat" *)

let load_lines filename =
  List.rev (Xlist.fold (Str.split (Str.regexp "\n") (File.load_file filename)) [] (fun l s ->
    if String.length s = 0 then s :: l else
    if String.sub s 0 1 = "#" then l else s :: l))

let fold_file filename n s f =
  let file = open_in filename in
  let r = ref s in
  try
    for i = 1 to n do
(*       if i mod 1000 = 0 then print_endline (string_of_int i);   *)
      let line = input_line file in
      r := f (!r) line
    done;
    close_in file;
    !r
  with
    End_of_file ->
      Printf.printf "file ended%!";
      close_in file;
      !r

let string_of_record_field_fmt name l =
  if l = [] then "" else
  name ^ " " ^ String.concat " " l ^ "\n"

let string_of_record_fmt r =
  r.rxn_id ^ " " ^ r.patent_id ^ "\n" ^
  r.reaction_smile ^ "\n" ^
  string_of_record_field_fmt "reactants" r.reactant_smiles ^
  string_of_record_field_fmt "products " r.product_smiles ^
(*  string_of_record_field_fmt "aux reactants" r.aux_reactant_smiles ^
  string_of_record_field_fmt "aux products " r.aux_product_smiles ^  *)
  string_of_record_field_fmt "solvents " r.solvent_smiles ^
(*   string_of_record_field_fmt "aux solvents " r.aux_solvent_smiles ^ *)
  "experimental yield " ^ r.experimental_yield ^ "  mintemp " ^ r.min_temp ^ "  max temp " ^ r.max_temp ^ "\n"

let rec get_list_stoi n = function
    x :: "" :: l ->
      if n = 0 then [],x :: "" :: l else
      let a,b = get_list_stoi (n-1) l in
      x :: a, b
  | l -> if n = 0 then [],l else failwith "get_list_stoi"

let rec get_list_stoi2 n = function
    x :: "" :: "" :: y :: "" :: "" :: l as ll ->
      if n = 0 then [],ll else
      let a,b = get_list_stoi2 (n-1) l in
      ((*map_stoi*) x,y) :: a, b
  | l -> if n = 0 then [],l else failwith "get_list_stoi2"

let rec get_list_rxn n = function
    x :: l ->
      if n = 0 then [],x :: l else
      let a,b = get_list_rxn (n-1) l in
      x :: a, b
  | [] -> if n = 0 then [],[] else failwith "get_list_rxn"

(*let map_stoi = function
  | "0.5" -> Half 1
  | "1" -> Aux 1
  | "1.0" -> Aux 1
  | "1.5" -> Half 3
  | "2" -> Aux 2
  | "2.0" -> Aux 2
  | "2.5" -> Half 5
  | "3" -> Aux 3
  | "3.0" -> Aux 3
  | "4" -> Aux 4
  | "4.0" -> Aux 4
  | "4.5" -> Half 9
  | "5" -> Aux 5
  | "5.0" -> Aux 5
  | "6" -> Aux 6
  | "6.0" -> Aux 6
  | "8" -> Aux 8
  | "8.0" -> Aux 8
  | "12.0" -> Aux 12
  | s -> failwith ("map_stoi: " ^ s)

let split_line_rxn_stoi n line =
  match Str.split (Str.regexp " ") line with
    rxn_id :: "" :: "" :: patent_id :: "" :: "" :: no_reactants :: "" :: "" :: no_products :: "" :: "" :: reaction_smile :: "" :: "" :: "" :: "" :: l ->
      let no_reactants = try int_of_string no_reactants with _ -> Printf.printf "%d: %s\n" n line; failwith "split_line_stoi 2" in
      let no_products = try int_of_string no_products with _ -> Printf.printf "%d: %s\n" n line; failwith "split_line_stoi 3" in
      let reactant_smiles,l = get_list_stoi no_reactants l in
      let product_smiles,l = get_list_stoi no_products l in
      (match l with
         "" :: no_missing_reactants :: "" :: "" :: no_missing_products :: "" :: "" :: l ->
           let no_missing_reactants = try int_of_string no_missing_reactants with _ -> Printf.printf "%d: %s\n" n line; failwith "split_line_stoi 5" in
           let no_missing_products = try int_of_string no_missing_products with _ -> Printf.printf "%d: %s\n" n line; failwith "split_line_stoi 6" in
           let missing_reactant_smiles,l = get_list_stoi2 no_missing_reactants l in
           let missing_product_smiles,l = get_list_stoi2 no_missing_products l in
           (match l with
              no_solvents :: "" :: l ->
                let no_solvents = try int_of_string no_solvents with _ -> Printf.printf "%d: %s\n" n line; failwith "split_line_stoi 7" in
                let solvents_smiles,l = get_list_stoi no_solvents l in
                (match l with
                   [""; experimental_yield; ""; ""; min_temp; ""; ""; max_temp] ->
                     {rxn_id = rxn_id; patent_id = patent_id;
                      reactant_smiles = reactant_smiles;
                      product_smiles = product_smiles;
                      solvent_smiles = solvents_smiles;
                      solvents = "";
                      missing_reactant_smiles = Xlist.map missing_reactant_smiles (fun (x,s) -> map_stoi x, s);
                      missing_product_smiles = Xlist.map missing_product_smiles (fun (x,s) -> map_stoi x, s);
                      aux_reactant_smiles = [];
                      aux_product_smiles = [];
                      aux_solvent_smiles = [];
                      reaction_smile = reaction_smile;
                      experimental_yield = experimental_yield;  experimental_yield_range = ""; experimental_yield_calc = "";
                      temp = ""; min_temp = min_temp; max_temp = max_temp;
                      bibref = ""}
                 | _ -> Printf.printf "%d: %s\n" n line; failwith "split_line_stoi 9")
(*                print_endline (String.concat ";" l);
                ()*)
            | _ -> Printf.printf "%d: %s\n" n line; failwith "split_line_stoi 8")
       | _ -> Printf.printf "%d: %s\n" n line; failwith "split_line_stoi 4")
  | _ -> Printf.printf "%d: %s\n" n line; failwith "split_line_stoi 1"*)

let split_line_rxn n line =
  match Str.split (Str.regexp "\t") line with
    rxn_id :: patent_id :: no_reactants :: no_products :: reaction_smile :: l ->
      let no_reactants = try int_of_string no_reactants with _ -> Printf.printf "%d: %s\n" n line; failwith "split_line_rxn 2" in
      let no_products = try int_of_string no_products with _ -> Printf.printf "%d: %s\n" n line; failwith "split_line_rxn 3" in
      let reactant_smiles,l = get_list_rxn no_reactants l in
      let product_smiles,l = get_list_rxn no_products l in
      (match l with
         no_solvents :: l ->
           let no_solvents = try int_of_string no_solvents with _ -> Printf.printf "%d: %s\n" n line; failwith "split_line_rxn 4" in
           let solvents_smiles,l = get_list_rxn no_solvents l in
           (match l with
              [experimental_yield; min_temp; max_temp] ->
                     {empty_record with
                      rxn_id = rxn_id; patent_id = patent_id;
                      reactant_smiles = reactant_smiles;
                      product_smiles = product_smiles;
                      solvent_smiles = solvents_smiles;
                      reaction_smile = reaction_smile;
                      experimental_yield = experimental_yield;
                      min_temp = min_temp; max_temp = max_temp}
            | _ -> Printf.printf "%d: %s\n" n line; failwith "split_line_stoi 6")
       | _ -> Printf.printf "%d: %s\n" n line; failwith "split_line_rxn 5")
  | _ -> Printf.printf "%d: %s\n" n line; failwith "split_line_rxn 1"

let fold_rxn s f =
  fold_file rxn_filename 1036446 s (fun s line ->
    let record = split_line_rxn 0 line in
    f s record)

let fold_rxn2 filename s f =
  fold_file filename 1036446 s (fun s line ->
    let record = split_line_rxn 0 line in
    f s record)

let fold_rxn_bounded n s f =
  fold_file rxn_filename n s (fun s line ->
    let record = split_line_rxn 0 line in
    f s record)

let load_rxn () =
  List.rev (fold_rxn [] (fun l r -> r :: l))

let load_rxn2 filename  =
  List.rev (fold_rxn2 filename [] (fun l r -> r :: l))

(*let fold_rxn_stoi s f =
  fold_file rxn_stoi_filename 259166 s (fun s line ->
    let record = split_line_rxn_stoi 0 line in
    f s record)*)

(*let rec process_paths name rev rev2 = function
    [] -> List.rev ((name,List.rev rev) :: rev2)
  | s :: l ->
      if String.length s = 0 then process_paths "" [] ((name,List.rev rev) :: rev2) l else
      if String.length s < 4 then failwith "process_paths" else
      if String.sub s 0 4 = "Path" then process_paths s rev rev2 l else
      process_paths name (s :: rev) rev2 l     *)

let load_paths () =
  let l = load_lines paths_filename in
(*   let l = process_paths "" [] [] (List.rev l) in
  Xlist.rev_map l (fun (name,l) ->
    name,*)
  List.rev (Xlist.rev_map l (fun s ->
    match Str.split (Str.regexp "\t") s with
      [(*path;*)id;smile] -> {empty_record with rxn_id = id; (*path = path;*) reaction_smile = smile}
    | _ -> failwith ("load_paths: " ^ s)))

let load_organic_syntheses () =
  let l = load_lines organic_syntheses_filename in
  List.rev (Xlist.rev_map l (fun s ->
    match Str.split (Str.regexp "\t") s with
      [id;smile] -> {empty_record with rxn_id = id; reaction_smile = smile}
    | _ -> failwith ("load_organic_syntheses: " ^ s)))

let load_organic_syntheses2 () =
  let l = load_lines organic_syntheses2_filename in
  List.rev (Xlist.rev_map l (fun s ->
    match Str.split (Str.regexp "\t") s with
      [id;smile;solvents;temp;y_range;y_avg;y_calc;bibref] ->  {empty_record with rxn_id = id;
                      solvents = solvents;
                      reaction_smile = smile;
                      experimental_yield = y_avg; experimental_yield_range = y_range; experimental_yield_calc = y_calc;
                      temp = temp;
                      bibref = bibref}
    | _ -> failwith ("load_organic_syntheses2: " ^ s)))

let load_przyklady_trudne () =
  let l = load_lines przyklady_trudne_filename in
  List.rev (Xlist.rev_map l (fun s ->
    match Str.split (Str.regexp "\t") s with
      [id;smile] -> {empty_record with rxn_id = id; reaction_smile = smile}
    | _ -> failwith ("load_przyklady_trudne: " ^ s)))

let load_przyklady_trudne2 () =
  let l = load_lines przyklady_trudne2_filename in
  List.rev (Xlist.rev_map l (fun s ->
    match Str.split (Str.regexp "\t") s with
      [id;smile] -> {empty_record with rxn_id = id; reaction_smile = smile}
    | _ -> failwith ("load_przyklady_trudne2: " ^ s)))

let load_losowe_patenty () =
  let l = load_lines losowe_patenty_filename in
  List.rev (Xlist.rev_map l (fun s ->
    match Str.split (Str.regexp "\t") s with
      [id;smile] -> {empty_record with rxn_id = id; reaction_smile = smile}
    | _ -> failwith ("load_losowe_patenty: " ^ s)))

let load_rt () =
  let l = File.load_tab rt_filename (function
    smile :: _ -> {empty_record with reaction_smile = smile}
  | l -> failwith ("load_rt: " ^ String.concat "\t" l)) in
  List.rev (snd (Xlist.fold (List.tl l) (1,[]) (fun (n,l) r ->
    n+1, {r with rxn_id = string_of_int n} :: l)))

let load_reactions filename =
  File.load_tab filename (function
    id :: smile :: _ -> {empty_record with rxn_id = id; reaction_smile = smile}
  | l -> failwith ("load_reactions: " ^ String.concat "\t" l))

let load_metacyc () =
  List.flatten (Xlist.map (Array.to_list (Sys.readdir metacyc_dir)) (fun filename ->
    load_reactions (metacyc_dir ^ filename)))

let string_of_is_correct = function
    true -> "YES"
  | false -> "NO"
    
let parse_is_correct = function
    "tak" -> true
  | "YES" -> true
  | "Tak" -> true
  | "nie" -> false
  | "NO (NO MAPPING)" -> false
  | "nie " -> false
  | "NO" -> false
  | "no" -> false
  | " nie" -> false
  | "Nie" -> false
  | "no matching found" -> false
  | s -> failwith ("parse_is_correct: '" ^ s ^ "'")
    
let parse_dataset = function
    "COMPLEX REACTIONS " -> "COMPLEX REACTIONS"
  | "trudne" -> "trudne"
  | "SIMPLE REACTIONS WITHOUT FULL STOICHIOMETRY" -> "SIMPLE REACTIONS WITHOUT FULL STOICHIOMETRY"
  | "SIMPLE REACTIONS WITH FULL STOICHIOMETRY" -> "SIMPLE REACTIONS WITH FULL STOICHIOMETRY"
  | "standardt" -> "standard"
  | "patenty" -> "patenty"
  | s -> failwith ("parse_dataset: '" ^ s ^ "'")
    
let load_wszystkie () =
  let l = load_lines wszystkie_filename in
  List.rev (Xlist.rev_map l (fun s ->
    match Str.split (Str.regexp ",") s with
      [id;smile;is_correct;dataset] -> {empty_record with rxn_id = id; reaction_smile = smile; path=parse_dataset dataset; is_correct=parse_is_correct is_correct}
    | _ -> failwith ("load_wszystkie: " ^ s)))

let create_molecules_list_fun set record =
  let reactant_smiles,product_smiles = Smiles.parse_smile_reaction record.reaction_smile in
  let set = Xlist.fold reactant_smiles set StringSet.add in
  let set = Xlist.fold product_smiles set StringSet.add in
  let set = Xlist.fold record.solvent_smiles set StringSet.add in
  set

let create_rxn_molecules_list filename =
  let molecules = fold_rxn StringSet.empty create_molecules_list_fun in
(*   let molecules = fold_rxn_stoi molecules create_molecules_list_fun in *)
  let molecules = StringSet.fold molecules StringSet.empty (fun set s ->
    let l = Str.split (Str.regexp "") s in
    if Smiles.is_error l then set else StringSet.add set s) in
  File.file_out filename (fun file ->
    StringSet.iter molecules (Printf.fprintf file "%s\n"))

let create_molecules_list records filename =
  let molecules = Xlist.fold records StringSet.empty (fun set r ->
    let reactant_smiles,product_smiles = Smiles.parse_smile_reaction r.reaction_smile in
    Xlist.fold (reactant_smiles @ product_smiles) set StringSet.add) in
  let molecules = Xlist.fold records molecules (fun set r ->
    Xlist.fold r.solvent_smiles set StringSet.add) in
  File.file_out filename (fun file ->
    StringSet.iter molecules (Printf.fprintf file "%s\n"))

(** Wytworzenie listy molekuł **)
(* let _ = create_rxn_molecules_list rxn_molecules_filename
let _ = create_molecules_list (load_paths ()) paths_molecules_filename
let _ = create_molecules_list (load_organic_syntheses ()) organic_syntheses_molecules_filename
let _ = create_molecules_list (load_organic_syntheses2 ()) organic_syntheses2_molecules_filename*)

let load_molecules = load_lines

(*let load_molecules_visualization filename =
  let l = load_lines filename in
  Xlist.fold l StringMap.empty (fun map s ->
    match Str.split (Str.regexp " ") s with
      smile :: smile2 :: smile3 :: l ->
        let tree = Smiles.parse_smile_molecule "" smile2 in
        let tree2 = Smiles.atom_tree_of_smile_tree_simple (Smiles.parse_smile_molecule "" smile3) in
        let angles = Array.of_list (Xlist.map l int_of_string) in
        StringMap.add map smile (tree,tree2,angles)
    | s -> failwith ("load_molecules_visualization: " ^ (String.concat " " s)))*)

(*let load_simple_reaction_ids () =
  StringSet.of_list (load_lines simple_reaction_ids_filename)*)
