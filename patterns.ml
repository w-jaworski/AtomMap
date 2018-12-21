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

type token =
    Text of string
  | Number of int
  | Paren of token list
  | Bracet of token list
  | SqBra of token list
  | LParen | RParen | LBracet | RBracet | LSqBra | RSqBra
  | Semic | Colon | Minus | Eq
  | P of pattern_atom_props

let make_tokens = function
    "0" -> Number 0
  | "1" -> Number 1
  | "2" -> Number 2
  | "3" -> Number 3
  | "4" -> Number 4
  | "5" -> Number 5
  | "6" -> Number 6
  | "7" -> Number 7
  | "8" -> Number 8
  | "9" -> Number 9
  | "(" -> LParen
  | ")" -> RParen
  | "{" -> LBracet
  | "}" -> RBracet
  | "[" -> LSqBra
  | "]" -> RSqBra
  | ";" -> Semic
  | ":" -> Colon
  | "-" -> Minus
  | "=" -> Eq
  | s -> Text s

let rec merge_text_numbers = function
    Text s :: Text t :: l -> Text(s ^ t) :: merge_text_numbers l
  | Number n :: Number m :: l -> Number((10*n)+m) :: merge_text_numbers l
  | t :: l -> t :: merge_text_numbers l
  | [] -> []

let rec find_brackets = function
    LParen :: l ->
        let found,l = find_rbracket RParen [] l in
        (Paren found) :: (find_brackets l)
  | LBracet :: l ->
        let found,l = find_rbracket RBracet [] l in
        (Bracet found) :: (find_brackets l)
  | LSqBra :: l ->
        let found,l = find_rbracket RSqBra [] l in
        (SqBra found) :: (find_brackets l)
  | s :: l -> s :: (find_brackets l)
  | [] -> []

and find_rbracket bracket rev = function
    LParen :: l ->
        let found,l = find_rbracket RParen [] l in
        find_rbracket bracket (Paren found :: rev) l
  | LBracet :: l ->
        let found,l = find_rbracket RBracet [] l in
        find_rbracket bracket (Bracet found :: rev) l
  | LSqBra :: l ->
        let found,l = find_rbracket RSqBra [] l in
        find_rbracket bracket (SqBra found :: rev) l
  | RParen :: l -> if bracket = RParen then List.rev rev, l else failwith "find_rbracket"
  | RBracet :: l -> if bracket = RBracet then List.rev rev, l else failwith "find_rbracket"
  | RSqBra :: l -> if bracket = RSqBra then List.rev rev, l else failwith "find_rbracket"
  | s :: l -> find_rbracket bracket (s :: rev) l
  | [] -> failwith "find_rbracket"

let rec parse_atom_names = function
    Text s :: Semic :: l -> s :: parse_atom_names l
  | [Text s] -> [s]
  | _ -> failwith "parse_atom_names"

let rec parse_atom_props s = function
    SqBra[Text name;Colon;Number id] :: l -> P{names=[name];pid=id} :: parse_atom_props s l
  | SqBra[Bracet names;Colon;Number id] :: l -> P{names=parse_atom_names names;pid=id} :: parse_atom_props s l
  | Paren p :: l -> Paren(parse_atom_props s p) :: parse_atom_props s l
  | t :: l -> t :: parse_atom_props s l
  | [] -> []

let rec parse_connections s = function
    Minus :: l  ->
      let atom,l = parse_atom s l in
      [Single,atom],l
  | Eq :: l  ->
      let atom,l = parse_atom s l in
      [Double,atom],l
  | Paren p :: l ->
      let conns,p = parse_connections s p in
      if p <> [] then failwith ("parse_connections: " ^ s) else
      let conns2,l = parse_connections s l in
      conns @ conns2,l
  | Number id :: l ->
      let conns,l = parse_connections s l in
      [Single,PLink id] @ conns,l
  | [] -> [],[]
  | _ -> failwith ("parse_connections: " ^ s)

and parse_atom s = function
    P p :: l ->
      let conns,l = parse_connections s l in
      PAtom(p,conns),l
  | Number id :: l -> PLink id,l
  | _ -> failwith ("parse_atom: " ^ s)


let parse_pattern s =
(*   let l = Xlist.map (Str.full_split (Str.regexp "]\|[:\|=\|-\|{\|}\|;\|(\|)") s) (function Str.Text s -> s | Str.Delim s -> s) in *)
  let l = Str.split (Str.regexp "") s in
  let l = Xlist.map l make_tokens in
  let l = merge_text_numbers l in
  let l = find_brackets l in
  let l = parse_atom_props s l in
  let t,l = parse_atom s l in
  if l <> [] then failwith ("parse_pattern: " ^ s) else
  t

let rec make_links_map links p = function
    _,PAtom(p,l) -> Xlist.fold l links (fun links a -> make_links_map links p a)
  | b,PLink t -> IntMap.add_inc links t (b,p,[]) (function (b2,p2,[]) -> b2,p2,[b,p] | _ -> failwith "make_links_map")

let add_edge graph from_node to_node bound =
  if (fst (graph.(from_node.pid))).pid = 0 then
    graph.(from_node.pid) <- from_node,[bound,to_node]
  else
    let _,l = graph.(from_node.pid) in
    graph.(from_node.pid) <- from_node,(bound,to_node) :: l

let rec make_pattern_graph_rec graph q = function
    b,PAtom(p,l) ->
      add_edge graph p q b;
      add_edge graph q p b;
      Xlist.iter l (fun a -> make_pattern_graph_rec graph p a)
  | _,PLink t -> ()

let match_bounds t = function
    Default,Single ->  Single
  | Single,Default -> Single
  | Default,x -> x (* FIXME: tu powinien byc blad *)
  | x,Default -> x
  | _ -> failwith ("match_bounds: " ^ string_of_int t)

let make_pattern_graph size l =
  let graph = Array.make size ({names=[]; pid=0},[]) in
  Xlist.iter l (function
      PAtom(p,l) ->
          graph.(p.pid) <- p,[];
          Xlist.iter l (fun a -> make_pattern_graph_rec graph p a);
          let links = Xlist.fold l IntMap.empty (fun links a -> make_links_map links p a) in
          IntMap.iter links (fun t -> function
              (b,p,[b2,q]) ->
                 let b = if b = b2 then b else match_bounds t (b,b2) in
                 add_edge graph p q b;
                 add_edge graph q p b
            | _ -> failwith "make_pattern_graph 3")
    | PLink t -> ()(*failwith ("make_atom_graph 1: " ^ t)*));
  graph

let rearrangement_patterns = [
  "Diels–Alder reaction",["[X:1]=[X:2]-[X:3]=[X:4]";"[X:5]=[X:6]"];
  "Sigmatropic",["[X:1]=[X:2]-[X:3]-[X:4]-[X:5]=[X:6]"];
  "Prins cyclization",["[C:1]=[C:2]-[C:3]-[C:4]-[O:5]";"[{C;H}:7]-[C:6](=[O:8])-[{C;H}:9]"];
  "Prins cyclization 2",["[C:2]=[C:3]-[C:5]1-[C:11]-[O:10]-[C:8]-[O:7]1"];
  "Wagner-Meerwein rearrangement",["[C:13]=[C:12]-[C:9]-[C:6]-[C:5]=[C:36]"];
  ]

let process_rearrangement_patterns l =
  Xlist.map l (fun (name,smarts) ->
    let _(*patterns*) = Xlist.map smarts (fun smart -> {smart=smart; smart_tree=parse_pattern smart}) in
    ())

(* let _ = process_rearrangement_patterns rearrangement_patterns *)
