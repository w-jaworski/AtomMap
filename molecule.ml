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

let rec get_link_paths path map = function
  | Atom(p,l) -> Xlist.fold l map (fun map (_,a) -> get_link_paths (p.id :: path) map a)
  | Link t -> StringMap.add_inc map t [List.rev path] (fun l -> (List.rev path) :: l)

let rec extract_cycles_rec t s = function
    n1 :: m1 :: p1, n2 :: m2 :: p2 ->
      if n1 <> n2 then failwith "extract_cycles_rec 1a" else
      if m1 = m2 then extract_cycles_rec t s (m1 :: p1,m2 :: p2) else
      [n1;m1;m2] @ p1 @ p2
  | [n1], n2 :: p2 ->
      if n1 <> n2 then failwith "extract_cycles_rec 1b" else
      [n1] @ p2
  | n1 :: p1, [n2] ->
      if n1 <> n2 then failwith "extract_cycles_rec 1c" else
      [n2] @ p1
  | p1,p2 ->
      Printf.printf "extract_cycles: bad paths [%s] [%s] to link %s in %s\n%!"
        (String.concat ";" (Xlist.map p1 string_of_int)) (String.concat ";" (Xlist.map p2 string_of_int)) t s;
      failwith "extract_cycles_rec 2"

let extract_cycles s paths =
  StringMap.fold paths [] (fun cycles t -> function
      [p1;p2] -> (extract_cycles_rec t s (p1,p2)) :: cycles
    | _ -> Printf.printf "extract_cycles: bad link %s in %s\n%!" t s; cycles)

let extract_labelled_cycles s paths =
  StringMap.mapi paths (fun t -> function
      [p1;p2] -> (extract_cycles_rec t s (p1,p2))
    | _ -> Printf.printf "extract_labelled_cycles: bad link %s in %s\n%!" t s; [])

let transitive_closure map =
  let f = ref true in
  let r = ref map in
  while !f do
    f := false;
    r := StringMap.fold (!r) (!r) (fun map k set ->
      let set2 = StringSet.fold set set (fun set2 v ->
        StringSet.union set2 (StringMap.find map v)) in
      if StringSet.size set2 > StringSet.size set then f := true;
      StringMap.add map k set2)
  done;
  !r

let cycle_union cycles =
  let map = StringMap.fold cycles IntMap.empty (fun map t l ->
    Xlist.fold l map (fun map n ->
      IntMap.add_inc map n [t] (fun l -> t :: l))) in
  let labels = IntMap.fold map StringMap.empty (fun labels _ l ->
    let set = Xlist.fold l StringSet.empty StringSet.add in
    Xlist.fold l labels (fun labels v ->
      StringMap.add_inc labels v set (fun set2 -> StringSet.union set set2))) in
  let labels = transitive_closure labels in
  let labels,_ = StringMap.fold labels ([],StringSet.empty) (fun (l,sum) v set ->
    if StringSet.mem sum v then l,sum else
    set :: l, StringSet.union sum set) in
  Xlist.map labels (fun set ->
    set, StringSet.fold set IntSet.empty (fun set v ->
      Xlist.fold (StringMap.find cycles v) set IntSet.add))

let get_molecule_core t =
  let paths = get_link_paths [] StringMap.empty t in
  let cycles = extract_labelled_cycles "" paths in
  let cycles = cycle_union cycles in
  IntSet.to_list (snd (Xlist.fold cycles (0,IntSet.empty) (fun (n,l) (links,nodes) ->
    if StringSet.size links > n then StringSet.size links, nodes else n,l)))

let rec stratify_molecule_rec cycles a n = function
    Atom(p,q) ->
      let n = if a.(p.id) = 0 then (
        let set = Xlist.fold cycles (IntSet.singleton p.id) (fun set (_,cycle) ->
          if IntSet.mem cycle p.id then cycle else set) in
        IntSet.iter set (fun i -> a.(i) <- n);
        n+1)
        else n in
      Xlist.iter q (fun (_,x) -> stratify_molecule_rec cycles a n x)
  | Link _ -> ()

let stratify_molecule size t =
  let a = Array.make size 0 in
  let paths = get_link_paths [] StringMap.empty t in
  let cycles = extract_labelled_cycles "" paths in
  let cycles = cycle_union cycles in
  stratify_molecule_rec cycles a 1 t;
  a



(****
open Types
open Xstd

let rec of_smile_rec m = function
  | b,SAtom(p,l) -> Path(Bond(m,b,p.id), Atom(p.id,p.name,p),Xlist.map l (of_smile_rec p.id))
  | b,SLink t -> Link(Bond(m,b,-1), t)
(*   | Default,SIon s -> Ion s *)
(*   | _,SIon _ -> failwith "of_smile_rec" *)

let of_smile a =
  match of_smile_rec (-1) (Default,a) with
    Path(_,node,paths) -> Root(node,paths)
  | _ -> failwith "of_smile"

let rec string_of_node = function
(*     Atom(_,s,_) -> s   *)
    Atom(n,s,_) -> s ^ string_of_int n
  | Complex(nodes,edges) ->
       "(" ^ String.concat "," (Xlist.map nodes string_of_node) ^ "|" ^
       String.concat "," (Xlist.map edges string_of_edge) ^ ")"

and string_of_path = function
    Path(edge,node,paths) ->
       "P(" ^ string_of_edge edge ^ "|" ^
       string_of_node node ^
       String.concat "" (Xlist.map paths string_of_path) ^ ")"
  | Link(edge,s) -> "L(" ^ string_of_edge edge ^ (if String.length s = 1 then s else "'" ^ s ^ "'") ^ ")"
  | Loop edge -> "O" ^ string_of_edge edge ^ ""
  | Leaf(edge,node,paths) ->
       "Le(" ^ string_of_edge edge ^ "|" ^
       string_of_node node ^
       String.concat "" (Xlist.map paths string_of_path) ^ ")"
  | Root(node,paths) -> string_of_node node ^ String.concat "" (Xlist.map paths string_of_path)
  | Ion s -> s

and string_of_edge = function
(*     Bond _ -> ""   *)
    Bond(m,_,n) -> Printf.sprintf "(%d-%d)" m n
  | Edge(l,e) -> "[" ^ String.concat "" (Xlist.map l (fun (e,a,paths) ->
      string_of_edge e ^ string_of_node a ^
      String.concat "" (Xlist.map paths string_of_path))) ^ string_of_edge e ^ "]"
  | Parallel l -> "{" ^ String.concat ";" (Xlist.map l string_of_edge) ^ "}"

(**********************************************************************************)
(** latex **)

let latex_of_bond = function
    Default -> ""
  | Single -> "-"
  | Slash -> "/"
  | Backslash -> "\\backslash "
  | Double -> "="
  | Triple -> "\\#"

let verbose_node = function
   Atom _ -> false
 | Complex _ -> true

let rec latex_of_node verbose = function
    Atom(n,s,_) -> if verbose then "\\text{{\\bf " ^ s ^ "}}_{" ^ string_of_int n ^ "}" else "\\text{{\\bf " ^ s ^ "}}"
  | Complex(nodes,edges) ->
       "{\\left(\\begin{array}{l}" ^ String.concat "\\\\\n" (Xlist.map nodes (latex_of_node true)) ^ "\\end{array}\\begin{array}{l}" ^
       String.concat "\\\\\n" (Xlist.map edges (latex_of_edge true true)) ^ "\\end{array}\\right)}"

and latex_of_path verbose = function
    Path(edge,node,paths) ->
       "\\text{{\\sc path}}& " ^ latex_of_edge verbose (verbose_node node) edge ^
       latex_of_node false node ^
       latex_of_paths (verbose_node node) paths
  | Link(edge,s) -> "\\text{{\\sc link}}& " ^ latex_of_edge verbose false edge ^ (if String.length s = 1 then s else "'" ^ s ^ "'")
  | Loop edge -> "\\text{{\\sc loop}}& " ^ latex_of_edge verbose verbose edge
  | Leaf(edge,node,paths) ->
       "\\text{{\\sc leaf}}& " ^ latex_of_edge verbose (verbose_node node) edge ^
       latex_of_node false node ^
       latex_of_paths (verbose_node node) paths
  | Root(node,paths) -> latex_of_node false node ^ latex_of_paths (verbose_node node) paths
  | Ion s -> "\\text{{\\sc ion}}& " ^ s

and latex_of_edge lv rv = function
    Bond(m,b,n) ->
        (match lv,rv with
          true,true -> Printf.sprintf "\\text{\\begin{scriptsize}%d\\end{scriptsize}}%s\\text{\\begin{scriptsize}%d\\end{scriptsize}}" m (latex_of_bond b) n
        | true,false -> Printf.sprintf "\\text{\\begin{scriptsize}%d\\end{scriptsize}%s}" m (latex_of_bond b)
        | false,true -> Printf.sprintf "%s\\text{\\begin{scriptsize}%d\\end{scriptsize}}" (latex_of_bond b) n
        | false,false -> latex_of_bond b)
  | Edge(l,e) ->
      let lv,s = Xlist.fold l (lv,"") (fun (lv,s) (e,a,paths) ->
        verbose_node a,
        s ^ latex_of_edge lv (verbose_node a) e ^ latex_of_node false a ^ latex_of_paths (verbose_node a) paths) in
      s ^ latex_of_edge lv rv e
  | Parallel l -> "{\\left\\{\\begin{array}{l}" ^ String.concat "\\\\\n" (Xlist.map l (latex_of_edge lv rv)) ^ "\\end{array}\\right\\}}"

and latex_of_paths verbose = function
    [] -> ""
  | paths -> "{\\left[\\begin{array}{ll}" ^ String.concat "\\\\\n" (Xlist.map paths (latex_of_path verbose)) ^ "\\end{array}\\right]}"

let latex_of t =
  "\\begin{flalign*}\n" ^ latex_of_path false t ^ "& &\\end{flalign*}\n"

let rec coloured_latex_of_node col verbose = function
    Atom(n,s,_) ->
      (if IntMap.mem col n then "\\textcolor{" ^ IntMap.find col n ^ "}{\\text{{\\bf " ^ s ^ "}}}" else "\\text{{\\bf " ^ s ^ "}}") ^
      if verbose then "_{" ^ string_of_int n ^ "}" else ""
  | Complex(nodes,edges) ->
       "{\\left(\\begin{array}{l}" ^ String.concat "\\\\\n" (Xlist.map nodes (coloured_latex_of_node col true)) ^ "\\end{array}\\begin{array}{l}" ^
       String.concat "\\\\\n" (Xlist.map edges (coloured_latex_of_edge col true true)) ^ "\\end{array}\\right)}"

and coloured_latex_of_path col verbose = function
    Path(edge,node,paths) ->
       "\\text{{\\sc path}}& " ^ coloured_latex_of_edge col verbose (verbose_node node) edge ^
       coloured_latex_of_node col false node ^
       coloured_latex_of_paths col (verbose_node node) paths
  | Link(edge,s) -> "\\text{{\\sc link}}& " ^ coloured_latex_of_edge col verbose false edge ^ (if String.length s = 1 then s else "'" ^ s ^ "'")
  | Loop edge -> "\\text{{\\sc loop}}& " ^ coloured_latex_of_edge col verbose verbose edge
  | Leaf(edge,node,paths) ->
       "\\text{{\\sc leaf}}& " ^ coloured_latex_of_edge col verbose (verbose_node node) edge ^
       coloured_latex_of_node col false node ^
       coloured_latex_of_paths col (verbose_node node) paths
  | Root(node,paths) -> coloured_latex_of_node col false node ^ coloured_latex_of_paths col (verbose_node node) paths
  | Ion s -> "\\text{{\\sc ion}}& " ^ s

and coloured_latex_of_edge col lv rv = function
    Bond(m,b,n) ->
        (match lv,rv with
          true,true -> Printf.sprintf "\\text{\\begin{scriptsize}%d\\end{scriptsize}}%s\\text{\\begin{scriptsize}%d\\end{scriptsize}}" m (latex_of_bond b) n
        | true,false -> Printf.sprintf "\\text{\\begin{scriptsize}%d\\end{scriptsize}%s}" m (latex_of_bond b)
        | false,true -> Printf.sprintf "%s\\text{\\begin{scriptsize}%d\\end{scriptsize}}" (latex_of_bond b) n
        | false,false -> latex_of_bond b)
  | Edge(l,e) ->
      let lv,s = Xlist.fold l (lv,"") (fun (lv,s) (e,a,paths) ->
        verbose_node a,
        s ^ coloured_latex_of_edge col lv (verbose_node a) e ^ coloured_latex_of_node col false a ^ coloured_latex_of_paths col (verbose_node a) paths) in
      s ^ coloured_latex_of_edge col lv rv e
  | Parallel l -> "{\\left\\{\\begin{array}{l}" ^ String.concat "\\\\\n" (Xlist.map l (coloured_latex_of_edge col lv rv)) ^ "\\end{array}\\right\\}}"

and coloured_latex_of_paths col verbose = function
    [] -> ""
  | paths -> "{\\left[\\begin{array}{ll}" ^ String.concat "\\\\\n" (Xlist.map paths (coloured_latex_of_path col verbose)) ^ "\\end{array}\\right]}"

let coloured_latex_of col t =
  "\\begin{flalign*}\n" ^ coloured_latex_of_path col false t ^ "& &\\end{flalign*}\n"

let rec edge_length = function
    Bond _ -> 0
  | Edge(l,e) ->
      Xlist.fold l 0 (fun n (e,a,_) -> n + 1 + edge_length e) + edge_length e
(*   | Parallel l -> "{\\left\\{\\begin{array}{l}" ^ String.concat "\\\\\n" (Xlist.map l (latex_chemfig_of_edge lv rv)) ^ "\\end{array}\\right\\}}" *)
  | _ -> failwith "edge_length: ni"

let rec latex_chemfig_of_node col verbose = function
    Atom(n,s,_) -> (*if verbose then "\\text{{\\bf " ^ s ^ "}}_{" ^ string_of_int n ^ "}" else "\\text{{\\bf " ^ s ^ "}}"*)
      (if IntMap.mem col n then "\\textcolor{" ^ IntMap.find col n ^ "}{\\text{{\\bf " ^ s ^ "}}}" else "\\text{{\\bf " ^ s ^ "}}") ^
      if verbose then "_{" ^ string_of_int n ^ "}" else ""
(*  | Complex(nodes,edges) ->
       "{\\left(\\begin{array}{l}" ^ String.concat "\\\\\n" (Xlist.map nodes (latex_chemfig_of_node col true)) ^ "\\end{array}\\begin{array}{l}" ^
       String.concat "\\\\\n" (Xlist.map edges (latex_chemfig_of_edge col true true)) ^ "\\end{array}\\right)}"*)
  | _ -> failwith "latex_chemfig_of_node: ni"

and latex_chemfig_of_path col angle verbose = function
(*    Path(edge,node,paths) ->
       "\\text{{\\sc path}}& " ^ latex_chemfig_of_edge col verbose (verbose_node node) edge ^
       latex_chemfig_of_node col false node ^
       latex_chemfig_of_paths col (verbose_node node) paths
  | Link(edge,s) -> "\\text{{\\sc link}}& " ^ latex_chemfig_of_edge col verbose false edge ^ (if String.length s = 1 then s else "'" ^ s ^ "'")*)
  | Loop edge -> Printf.sprintf "*%d(%s)" (edge_length edge + 1) (latex_chemfig_of_edge col "" verbose verbose edge)
  | Leaf(edge,node,paths) ->
       latex_chemfig_of_edge col angle verbose (verbose_node node) edge ^
       latex_chemfig_of_node col false node ^
       latex_chemfig_of_paths col 0 (verbose_node node) paths
  | Root(node,paths) -> latex_chemfig_of_node col false node ^ latex_chemfig_of_paths col 0 (verbose_node node) paths
(*  | Ion s -> "\\text{{\\sc ion}}& " ^ s*)
  | _ -> failwith "latex_chemfig_of_path: ni"

and latex_chemfig_of_edge col angle lv rv = function
    Bond(m,b,n) ->
        (Smiles.latex_chemfig_of_bond b) ^ angle
(*        (match lv,rv with
          true,true -> Printf.sprintf "\\text{\\begin{scriptsize}%d\\end{scriptsize}}%s\\text{\\begin{scriptsize}%d\\end{scriptsize}}" m (latex_chemfig_of_bond b) n
        | true,false -> Printf.sprintf "\\text{\\begin{scriptsize}%d\\end{scriptsize}%s}" m (latex_chemfig_of_bond b)
        | false,true -> Printf.sprintf "%s\\text{\\begin{scriptsize}%d\\end{scriptsize}}" (latex_chemfig_of_bond b) n
        | false,false -> latex_chemfig_of_bond b)*)
  | Edge(l,e) ->
      let lv,s,angle = Xlist.fold l (lv,"",angle) (fun (lv,s,angle) (e,a,paths) ->
        verbose_node a,
        s ^ latex_chemfig_of_edge col angle lv (verbose_node a) e ^ latex_chemfig_of_node col false a ^ latex_chemfig_of_paths col 1 (verbose_node a) paths,
        (*""*)angle) in
      s ^ latex_chemfig_of_edge col angle lv rv e
(*   | Parallel l -> "{\\left\\{\\begin{array}{l}" ^ String.concat "\\\\\n" (Xlist.map l (latex_chemfig_of_edge col lv rv)) ^ "\\end{array}\\right\\}}" *)
  | _ -> failwith "latex_chemfig_of_edge: ni"

and latex_chemfig_of_paths col x verbose = function
    [] -> ""
  | paths ->
      let n = Xlist.size paths + 3 in
      fst (Xlist.fold paths ("",x) (fun (s,i) path ->
        s ^ (if i<n-1+x then "(" else "") ^
        latex_chemfig_of_path col (Printf.sprintf "[::%d]" (360*i/(n+x))) verbose path ^
        (if i<n-1+x then ")" else ""),i+1))

let latex_chemfig_of col t =
  "\\chemfig{" ^ latex_chemfig_of_path col "" false t ^ "}\\\\"


(*let rec latex_chemfig_of_symbol d = function
    b,Aroma(s,t,a) -> latex_chemfig_of_symbol d (b,Atom(s ^ "^{*}",t,a))
  | Root,Atom(s,t,[]) -> s
  | Root,Atom(s,t,[x]) -> s ^ latex_chemfig_of_symbol "" x
  | Root,Atom(s,t,[x;y]) -> s ^ latex_chemfig_of_symbol "[::90]" x ^ latex_chemfig_of_symbol "[::0]" y
  | Root,Atom(s,t,[x;y;z]) -> s ^ latex_chemfig_of_symbol "[::90]" x ^ latex_chemfig_of_symbol "[::-90]" y ^ latex_chemfig_of_symbol "[::0]" z
  | Root,Atom(s,t,[h;x;y;z]) -> s ^ latex_chemfig_of_symbol "[::180]" h ^ latex_chemfig_of_symbol "[::90]" x ^ latex_chemfig_of_symbol "[::-90]" y ^ latex_chemfig_of_symbol "[::0]" z
  | b,Atom(s,t,[]) -> "(" ^ latex_chemfig_of_bond b ^ d ^ s ^ ")"
  | b,Atom(s,t,[x]) -> "(" ^ latex_chemfig_of_bond b ^ d ^ s ^ latex_chemfig_of_symbol "[::0]" x ^ ")"
  | b,Atom(s,t,[x;y]) -> "(" ^ latex_chemfig_of_bond b ^ d ^ s ^ latex_chemfig_of_symbol "[::90]" x ^ latex_chemfig_of_symbol "[::0]" y ^ ")"
  | b,Atom(s,t,[x;y;z]) -> "(" ^ latex_chemfig_of_bond b ^ d ^ s ^ latex_chemfig_of_symbol "[::90]" x ^ latex_chemfig_of_symbol "[::-90]" y ^ latex_chemfig_of_symbol "[::0]" z ^ ")"
  | b,Number s -> "(" ^ latex_chemfig_of_bond b ^ d ^ s ^ ")"
  | Root,AromaRing(Atom(s,t,[]) :: l) ->
      s ^ "**" ^ (string_of_int (1 + Xlist.size l)) ^
      "(-" ^ String.concat "-" (Xlist.map l (fun x -> latex_chemfig_of_symbol "" (Root,x))) ^ "-)"
  | b,AromaRing(Atom(s,t,[]) :: l) ->
      "(" ^ latex_chemfig_of_bond b ^ d ^ s ^ "**" ^ (string_of_int (1 + Xlist.size l)) ^
      "(-" ^ String.concat "-" (Xlist.map l (fun x -> latex_chemfig_of_symbol "" (Root,x))) ^ "-))"
  | _ -> failwith "latex_chemfig_of_symbol"*)

(**********************************************************************************)
(** ring recognition **)

let rec count_links_path set = function
    Path(_,_,paths) -> Xlist.fold paths set count_links_path
  | Link(_,s) -> StringSet.add set s
  | Root(_,paths) -> Xlist.fold paths set count_links_path
  | _ -> set

let count_links t =
  let set = count_links_path StringSet.empty t in
  StringSet.size set

let rec split_leaf_paths (a,b) = function
  | Leaf(edge, node, paths) :: l -> split_leaf_paths (Leaf(edge, node, paths) :: a,b) l
  | Root _ :: _ -> failwith "split_leaf_paths"
  | t :: l -> split_leaf_paths (a,t :: b) l
  | [] -> a,b

let rec split_trivial_paths (a,b) = function
  | Loop edge :: l -> split_trivial_paths (Loop edge :: a,b) l
  | Leaf(edge, node, paths) :: l -> split_trivial_paths (Leaf(edge, node, paths) :: a,b) l
  | Ion s :: l -> split_trivial_paths (Ion s :: a,b) l
  | Root _ :: _ -> failwith "split_trivial_paths"
  | t :: l -> split_trivial_paths (a,t :: b) l
  | [] -> a,b

let remove_trivial_paths l =
  snd (split_trivial_paths ([],[]) l)

let rec revert_edge = function
    Bond(m,b,n) -> Bond(n,b,m)
  | Edge(l,e) -> let l,e = revert_edge_list e (List.rev l) in Edge(l,e)
  | Parallel l -> Parallel(Xlist.map l revert_edge)

and revert_edge_list e_last = function
    [] -> [],revert_edge e_last
  | (e,a,p) :: l ->
      let l,e = revert_edge_list e l in
      (revert_edge e_last,a,p) :: l, e



let merge_edges a p = function
    Edge(l1,e1),Edge(l2,e2) -> Edge(l1 @ [e1,a,p] @ l2,e2)
  | e1,Edge(l2,e2) -> Edge([e1,a,p] @ l2,e2)
  | Edge(l1,e1),e2 -> Edge(l1 @ [e1,a,p],e2)
  | e1,e2 -> Edge([e1,a,p],e2)

let rec merge_edges_link = function
    Edge(l1,e1),Edge((e2,a,p) :: l2,e3) -> Edge(l1 @ [merge_edges_link (e1,e2),a,p] @ l2,e3)
  | Bond(-2,Default,-2),Bond(-2,Default,-2) -> failwith "merge_edges_link 6"
  | Bond(-2,Default,-2),e2 -> e2
  | e1,Bond(-2,Default,-2) -> e1
  | _, Edge([],_) -> failwith "merge_edges_link 4"
  | Edge([],_), _ -> failwith "merge_edges_link 5"
  | e1,Edge((e2,a,p) :: l2,e3) -> Edge([merge_edges_link (e1,e2),a,p] @ l2,e3)
  | Edge(l1,e1),e2 -> Edge(l1,merge_edges_link (e1,e2))
  | Bond(m,b1,_),Bond(_,b2,n) -> if b1 = b2 then Bond(m,b1,n) else failwith "merge_edges_link: bond"
  | Parallel _,_ -> failwith "merge_edges_link 1"
  | _, Parallel _ -> failwith "merge_edges_link 2"

let rec get_first_bond = function
    Edge((e,_,_) :: _,_) -> get_first_bond e
  | Bond(m,b,_) -> Bond(m,b,-1)
  | Parallel _ -> Bond(-2,Default,-2)
  | _ -> failwith "get_first_bond"

let rec merge_paths = function
    Path(edge,node,paths) ->
      let tpaths,ntpaths = split_trivial_paths ([],[]) paths in
      (match ntpaths with
        [Path(edge2,a,paths2)] -> merge_paths (Path(merge_edges node tpaths (edge,edge2),a,paths2))
      | [Link(edge2,s2)] -> Link(merge_edges node tpaths (edge,edge2),s2)
      | _ -> Path(edge,node,Xlist.map paths merge_paths))
  | Root(node,paths) -> Root(node,Xlist.map paths merge_paths)
  | path -> path

let rec mark_leaf_paths = function
    Path(edge,node,paths) ->
      let paths = Xlist.map paths mark_leaf_paths in
      let tpaths,ntpaths = split_trivial_paths ([],[]) paths in
      (match ntpaths with
        [] -> Leaf(edge,node,paths)
      | _ -> Path(edge,node,paths))
  | Root(node,paths) -> Root(node,Xlist.map paths mark_leaf_paths)
  | path -> path

let rec merge_leaf_paths = function
    Path(edge,node,paths) -> Path(edge,node,Xlist.map paths merge_leaf_paths)
  | Leaf(edge,node,paths) ->
      let lpaths,nlpaths = split_leaf_paths ([],[]) paths in
      (match lpaths with
        [] -> Leaf(edge,node,paths)
      | [Leaf(edge2,node2,paths2)] -> merge_leaf_paths (Leaf(merge_edges node nlpaths (edge,edge2),node2,paths2))
      | _ -> Leaf(edge,node,Xlist.map paths merge_leaf_paths))
  | Root(node,paths) -> Root(node,Xlist.map paths merge_leaf_paths)
  | path -> path

(*let rec move_paths map = function
    Atom(n,s,p,paths) ->
       let map,paths = Xlist.fold paths (map,[]) (fun (map,paths) path ->
         let map,path = move_paths_path map path in
         map, path :: paths) in
       map,Atom (n,s,p,paths)
  | Complex(nodes,edges,paths) ->
       let map,paths = Xlist.fold paths (map,[]) (fun (map,paths) path ->
         let map,path = move_paths_path map path in
         map, path :: paths) in
       map,Complex(nodes,edges,paths)

and move_paths_path map = function
    Path(edge,a) -> let map,a = move_paths map a in map,Path(edge,a)
  | Link(edge,t) ->
      if StringMap.mem map t then
        StringMap.remove map t, Link(merge_edges_link (edge,revert_edge (StringMap.find map t)), t)
      else StringMap.add map t edge, Link(get_first_bond edge, t)
  | Loop edge -> map, Loop edge
  | Leaf(edge, node) -> map, Leaf(edge, node)*)

let rec mark_loops_paths paths =
  let paths = Xlist.map paths mark_loops in
  let links,others = Xlist.fold paths (StringMap.empty,[]) (fun (links,others) -> function
           Link(edge,t) -> StringMap.add_inc links t [edge,t] (fun l -> (edge,t) :: l), others
         | a -> links, a :: others) in
  StringMap.fold links others (fun paths _ -> function
           [edge,t] -> Link(edge,t) :: paths
         | [edge1,_;edge2,_] -> Loop(merge_edges_link (edge1,revert_edge edge2)) :: paths
         | _ -> failwith "mark_loops_paths")

and mark_loops = function
     Path(edge,node,paths) -> Path(edge,node,mark_loops_paths paths)
   | Root(node,paths) -> Root(node,mark_loops_paths paths)
   | a -> a

let arrange_connectives (links,others) = function
    Link(edge,t) -> StringMap.add links t (Link(edge,t)), others
  | t -> links, t :: others

let get_edge = function
    Link(edge,_) -> edge
  | t -> failwith ("get_edge: " ^ string_of_path t)

let mark_parallels_paths2 links link_ids paths edge node paths2 =
  let links2,others2 = Xlist.fold paths2 (StringMap.empty,[]) arrange_connectives in
  let link_ids2 = StringMap.fold links2 StringSet.empty (fun set k _ -> StringSet.add set k) in
  let c = StringSet.intersection link_ids link_ids2 in
  if StringSet.size c = 0 then Path(edge,node,paths2) :: paths, links else
  let new_edges,links,links2 = StringSet.fold c ([],links,links2) (fun (new_edges,links,links2) id ->
        merge_edges_link (get_edge (StringMap.find links id),revert_edge (get_edge (StringMap.find links2 id))) :: new_edges,
        StringMap.remove links id,
        StringMap.remove links2 id) in
  Path(Parallel(edge :: new_edges),node,StringMap.fold links2 others2 (fun paths2 _ t -> t :: paths2)) :: paths, links

let rec mark_parallels_paths paths =
  let links,others = Xlist.fold paths (StringMap.empty,[]) arrange_connectives in
  let link_ids = StringMap.fold links StringSet.empty (fun set k _ -> StringSet.add set k) in
  let paths,links = Xlist.fold others ([],links) (fun (paths,links) -> function
           Path(edge,node,paths2) -> mark_parallels_paths2 links link_ids paths edge node paths2
         | t -> t :: paths, links) in
  let paths = StringMap.fold links paths (fun paths _ t -> t :: paths) in
  Xlist.map paths mark_parallels

and mark_parallels = function
  | Path(edge,node,paths) -> Path(edge,node,mark_parallels_paths paths)
  | Root(node,paths) -> Root(node,mark_parallels_paths paths)
  | t -> t

let rec mark_distant_parallels_paths2 links paths link_ids edge node paths2 =
  let links2,others2 = Xlist.fold paths2 (StringMap.empty,[]) arrange_connectives in
  let link_ids2 = StringMap.fold links2 StringSet.empty (fun set k _ -> StringSet.add set k) in
  let c = StringSet.intersection link_ids link_ids2 in
  if StringSet.size c < 2 then
    let paths,links,paths2 = Xlist.fold paths2 (paths,links,[]) (fun (paths,links,paths2) a ->
      let paths,links,a = mark_distant_parallels2_path links paths link_ids a in
      paths,links,a :: paths2) in
    paths, links, Path(edge,node,paths2)
  else
    let new_edges,links,links2 = StringSet.fold c ([],links,links2) (fun (new_edges,links,links2) id ->
      merge_edges_link (get_edge (StringMap.find links id),revert_edge (get_edge (StringMap.find links2 id))) :: new_edges,
      StringMap.remove links id,
      StringMap.remove links2 id) in
    let x = StringSet.min_elt c in
    let a = StringMap.fold links2 (Link(get_first_bond (revert_edge edge),x) :: others2) (fun paths2 _ t -> t :: paths2) in
    let paths,links,a = mark_distant_parallels_paths2 links paths link_ids (Parallel new_edges) node a in
    a :: paths, links, Link(edge,x)

and mark_distant_parallels2_path links paths link_ids = function
  | Path(edge,node,paths2) ->
     mark_distant_parallels_paths2 links paths link_ids edge node paths2
  | t -> paths, links, t

let rec mark_distant_parallels_paths paths =
  let links,others = Xlist.fold paths (StringMap.empty,[]) arrange_connectives in
  let link_ids = StringMap.fold links StringSet.empty (fun set k _ -> StringSet.add set k) in
  let paths,links = Xlist.fold others ([],links) (fun (paths,links) t ->
    let paths,links,t = mark_distant_parallels2_path links paths link_ids t in
    t :: paths, links) in
  let paths = StringMap.fold links paths (fun paths _ t -> t :: paths) in
  Xlist.map paths mark_distant_parallels

and mark_distant_parallels = function
  | Path(edge,node,paths) -> Path(edge,node,mark_distant_parallels_paths paths)
  | Root(node,paths) -> Root(node,mark_distant_parallels_paths paths)
  | t -> t

let rec mark_triangles3_path node links link_ids others edge node2 rev2 = function
  | Path(edge2,node3,paths3) :: paths2 ->
       (try
         let links3,others3 = Xlist.fold paths3 (StringMap.empty,[]) arrange_connectives in
         let link_ids3 = StringMap.fold links3 StringSet.empty (fun set k _ -> StringSet.add set k) in
         let c = StringSet.intersection link_ids link_ids3 in
         if StringSet.size c = 0 then raise Not_found else
         let id = StringSet.min_elt c in
         let paths = StringMap.fold (StringMap.remove links id) others (fun paths _ t -> t :: paths) in
         let paths3 = StringMap.fold (StringMap.remove links3 id) others3 (fun paths3 _ t -> t :: paths3) in
         Complex([node;node2;node3],
           [edge; edge2; merge_edges_link (get_edge (StringMap.find links3 id),revert_edge (get_edge (StringMap.find links id)))]),
           paths @ rev2 @ paths2 @ paths3
       with Not_found -> mark_triangles3_path node links link_ids others edge node2 (Path(edge2,node3,paths3) :: rev2) paths2)
  | t :: paths2 -> mark_triangles3_path node links link_ids others edge node2 (t :: rev2) paths2
  | [] -> raise Not_found

let rec mark_triangles2_path node links link_ids rev = function
  | Path(edge,node2,paths2) :: others ->
       (try mark_triangles3_path node links link_ids (rev @ others) edge node2 [] paths2
       with Not_found -> mark_triangles2_path node links link_ids (Path(edge,node2,paths2) :: rev) others)
  | t :: others -> mark_triangles2_path node links link_ids (t :: rev) others
  | [] -> raise Not_found

let rec mark_triangles_node node paths =
  let links,others = Xlist.fold paths (StringMap.empty,[]) arrange_connectives in
  let link_ids = StringMap.fold links StringSet.empty (fun set k _ -> StringSet.add set k) in
  (try mark_triangles2_path node links link_ids [] others with
         Not_found -> node, mark_triangles_path paths)

and mark_triangles_path = function
  | Path(edge,node,paths) :: path ->
      (try
         let node,paths = mark_triangles_node node paths in
         Path(edge,node,paths) :: path
       with Not_found -> Path(edge,node,paths) :: (mark_triangles_path path))
  | t :: path -> t :: mark_triangles_path path
  | [] -> raise Not_found

let mark_triangles = function
    Root(node,paths) ->
      let node,paths = mark_triangles_node node paths in
      Root(node,paths)
  | _ -> failwith "mark_triangles"

let rec process_loops t =
  let c = count_links t in
(*  let deg = find_max_degree t in
(*     Printf.printf "l1 %s\n%!" (string_of_path t);   *)
   let t = if deg > 2 then try List.hd (change_root deg t) with _ -> t else t in *)
(*     Printf.printf "l1a %s\n%!" (string_of_path t);   *)
  let t = mark_leaf_paths t in
  let t = merge_leaf_paths t in
  let t = merge_paths t in
(*   let _,t = move_paths StringMap.empty t in *)
(*    Printf.printf "l2 %s\n%!" (string_of_path t);   *)
  let t = mark_loops t in
(*    Printf.printf "l3 %s\n%!" (string_of_path t);   *)
  if c = count_links t then t else process_loops t

let rec process_parallels t =
(*      Printf.printf "p1 %s\n%!" (string_of_path t);     *)
  let c = count_links t in
(*      Printf.printf "p1a %s\n%!" (string_of_path t);     *)
  let t = process_loops t in
(*      Printf.printf "p2 %s\n%!" (string_of_path t);     *)
(*  let deg = find_max_degree t in
   let t = if deg > 2 then try List.hd (change_root deg t) with _ -> t else t in *)
  let t = mark_leaf_paths t in
  let t = merge_leaf_paths t in
  let t = merge_paths t in
(*   let _,t = move_paths StringMap.empty t in *)
(*      Printf.printf "p3 %s\n%!" (string_of_path t);     *)
  let t = mark_parallels t in
(*      Printf.printf "p4 %s\n%!" (string_of_path t);     *)
  let t = mark_distant_parallels t in
(*      Printf.printf "p5 %s\n%!" (string_of_path t);     *)
  if c = count_links t then t else process_parallels t

let rec process_triangles t =
(*      Printf.printf "t1 %s\n%!" (string_of_path t);     *)
  let c = count_links t in
  let t = process_parallels t in
(*      Printf.printf "t2 %s\n%!" (string_of_path t);     *)
(*  let deg = find_max_degree t in
   let t = if deg > 2 then try List.hd (change_root deg t) with _ -> t else t in *)
  let t = mark_leaf_paths t in
  let t = merge_leaf_paths t in
  let t = merge_paths t in
(*   let _,t = move_paths StringMap.empty t in *)
(*      Printf.printf "t3 %s\n%!" (string_of_path t);    *)
  let t = try mark_triangles t with Not_found -> t in
(*      Printf.printf "t4 %s\n%!" (string_of_path t);     *)
(*   let t = mark_distant_triangles t in    *)
  if c = count_links t then t else process_triangles t

(**********************************************************************************)
(** pattern creation **)

let rec mark_coloured_node col = function
    Atom(n,s,p) -> if IntMap.mem col n then Atom(n,"X",p) else Atom(n,s,p)
  | Complex(nodes,edges) -> Complex(Xlist.map nodes (mark_coloured_node col),Xlist.map edges (mark_coloured_edge col))

and mark_coloured_path col = function
    Path(edge,node,paths) -> Path(mark_coloured_edge col edge,mark_coloured_node col node,Xlist.map paths (mark_coloured_path col))
  | Link(edge,s) -> Link(mark_coloured_edge col edge,s)
  | Loop edge -> Loop(mark_coloured_edge col edge)
  | Leaf(edge,node,paths) -> Leaf(mark_coloured_edge col edge,mark_coloured_node col node,Xlist.map paths (mark_coloured_path col))
  | Root(node,paths) -> Root(mark_coloured_node col node,Xlist.map paths (mark_coloured_path col))
  | Ion s -> Ion s

and mark_coloured_edge col = function
    Bond(m,b,n) -> Bond(m,b,n)
  | Edge(l,e) -> Edge(Xlist.map l (fun (e,a,paths) ->
     mark_coloured_edge col e, mark_coloured_node col a, Xlist.map paths (mark_coloured_path col)),
     mark_coloured_edge col e)
  | Parallel l -> Parallel(Xlist.map l (mark_coloured_edge col))

let rec reduce_coloured_node = function
    Atom(n,s,p) -> Atom(n,s,p)
  | Complex(nodes,edges) -> Complex(Xlist.map nodes reduce_coloured_node,Xlist.map edges reduce_coloured_edge)

and reduce_coloured_path = function
    Path(edge,node,paths) -> Path(reduce_coloured_edge edge,reduce_coloured_node node,reduce_coloured_paths paths)
  | Link(edge,s) -> Link(reduce_coloured_edge edge,s)
  | Loop edge ->
       (match reduce_coloured_edge edge with
         Bond(_,_,_) -> Ion ""
       | edge -> Loop edge)
  | Leaf(edge,node,paths) ->
       (match reduce_coloured_edge edge,reduce_coloured_node node,reduce_coloured_paths paths with
         Bond(_,_,_),Atom(_,"X",_),[] -> Ion ""
       | edge,node,paths -> Leaf(edge,node,paths))
  | Root(node,paths) -> Root(reduce_coloured_node node,reduce_coloured_paths paths)
  | Ion s -> Ion s

and reduce_coloured_paths = function
    path :: paths ->
      (match reduce_coloured_path path with
        Ion "" -> reduce_coloured_paths paths
      | path -> path :: reduce_coloured_paths paths)
  | [] -> []

and reduce_coloured_edge = function
    Bond(m,b,n) -> Bond(m,b,n)
  | Edge(l,e) ->
      let l = reduce_coloured_edge_list (Xlist.map l (fun (e,a,paths) ->
        reduce_coloured_edge e, reduce_coloured_node a, reduce_coloured_paths paths)) in
      (match l,reduce_coloured_edge e with
        [Bond(m,b,_),Atom(_,"X",_),[]],Bond(_,_,n) -> Bond(m,b,n)
      | l,e -> Edge(l,e))
  | Parallel l -> Parallel(Xlist.map l (reduce_coloured_edge))

and reduce_coloured_edge_list = function
  | (Bond(m,b,_),Atom(_,"X",_),[]) :: (Bond(_,_,n),a,paths) :: l -> reduce_coloured_edge_list ((Bond(m,b,n),a,paths) :: l)
  | x :: y :: l -> x :: (reduce_coloured_edge_list (y :: l))
  | l -> l

let reduce_coloured col t =
  let t = mark_coloured_path col t in
  reduce_coloured_path t

****)
