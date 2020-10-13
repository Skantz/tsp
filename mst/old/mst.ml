#load "str.cma";;

module IntSet = Set.Make( 
  struct
    let compare = Pervasives.compare
    type t = int
  end)

let query_node_local_cycle n max_n node_lookup = 
  let () = assert (max_n > 2) in
  let x = (Hashtbl.find_opt node_lookup n) in
  match x with
  | None -> let () = Printf.printf "%d\n" n in
    let line = [n - 1 mod max_n; 1; n + 1 mod max_n; 1] in
    let () = Hashtbl.add node_lookup n (Some line) in
    line
  | Some x -> Option.get x
;;

let query_node n node_lookup =
  (*let n : int = Random.int max_n in*)
  let x = (Hashtbl.find_opt node_lookup n) in
  match x with
  | None -> let () = Printf.printf "%d\n" n in
    let line = read_line () in
    let line = Str.split (Str.regexp_string " ") line |> List.map int_of_string in
    let line = List.tl line in
    let () = Hashtbl.add node_lookup n (Some line) in
    line
  | Some x -> Option.get x
;;

let bfs start_node max_n node_lookup =
  let q = Queue.create () in
  let bound = int_of_float (1.0 /. Random.float 1.0) in
  let () = Queue.add start_node q in
  let i = ref 0 in
  let continue = ref true in
  let visited = ref (IntSet.empty) in
  while !i < bound && !continue do 
    let n = Queue.take q in
    (*let neighbors = query_node n node_lookup in*)
    let neighbors = query_node_local_cycle n max_n node_lookup in
    let neighbors = List.filter (fun x -> (IntSet.mem x !visited))  neighbors in
    if (neighbors = []) then continue := false; 
    List.iter (fun x -> Queue.add x q)         neighbors;
    List.iter (fun x -> (IntSet.add x !visited )) neighbors;
    i := (!i + 1);
  done;
  (!i = bound)

let est_mst n eps max_w node_lookup =
  let res = Array.init max_w (fun x -> 0) in 
  let bound = 1 + int_of_float ((float_of_int max_w) /. (eps ** 2.0)) in
  let i = ref 0 in 
  while !i < bound do
    let start_node = Random.int (n - 1) in
    let w = ref max_w in
    let bool_to_int t = if t then 1 else 0 in
    while !w > 0 do
      res.(!w - 1) <- res.(!w - 1) + bool_to_int (bfs start_node n node_lookup);
      w := !w - 1;
      ()
    done;
    ();
    i := !i + 1
  done;
  Array.map (fun x -> (x * n / bound)) res
;;
(*
let print_random_number max_n =
  let n = Random.int max_n in
  Printf.printf "%d\n" n
*)

Random.self_init ()

let n = read_int () 
let eps = read_float () -. 1.0 
let max_w = read_int ()

let node_lookup = Hashtbl.create 1000

(*print_random_number 15
let main () = query_node 10 node_lookup ()*)
let ans = est_mst n eps max_w node_lookup
let main = Array.iter (Printf.printf "%d ") ans