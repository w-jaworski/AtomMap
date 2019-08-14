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

let load_reactions filename =
  File.load_tab filename (function
    id :: smile :: _ -> {empty_record with rxn_id = id; reaction_smile = smile}
  | l -> failwith ("load_reactions: " ^ String.concat "\t" l))

let string_of_is_correct = function
    true -> "YES"
  | false -> "NO"
    
