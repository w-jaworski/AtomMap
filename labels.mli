(*
 *  AtomMap: maps atoms in chemical reactions
 *  Copyright (C) 2015-2017 Wojciech Jaworski <wjaworski atSPAMfree mimuw dot edu dot pl> 
 *  Copyright (C) 2015-2017 Institute of Informatics, University of Warsaw                

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


open Xstd

type t

(* val empty : t *)
val make : int -> t
val find : t -> int -> int
val get : t -> int -> int
(*val set : t -> int -> unit
val set_if_new : t -> int -> unit
val next : t -> unit*)
val get_next : t -> int
val mem : t -> int -> bool
val fold : t -> 'b -> ('b -> int -> 'b) -> 'b
val foldi : t -> 'b -> ('b -> int -> int -> 'b) -> 'b
val add_list_of_list_pairs : t -> int Pair.t -> t
(* val make_list_of_lists : int -> int list list -> t *)
val add_list_of_lists : t -> int list list -> t
val map : t -> int IntMap.t -> int -> t
val remove_list : t -> int list -> t
val add_invisible_list : t -> int list -> t
