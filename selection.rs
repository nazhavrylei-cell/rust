// Селекція рослин і тварин — консольна гра на Rust (без зовнішніх залежностей).
//
// Збірка і запуск:   rustc -O selection.rs && ./selection
//
// Модель генетики (за Менделем):
//   * кожен ген має два алелі: велика літера = домінантний, мала = рецесивний;
//   * геном записується парами літер: AaBbccDd (AA, Aa, aa);
//   * потомок від кожного батька отримує по одному випадковому алелю гена.
//
// Програма шукає всі генотипи, які можна отримати схрещуваннями з ваших особин,
// вибирає найкращий (за вагою генів) і показує найкоротший план схрещувань
// із шансами на кожному кроці.

use std::cmp::Ordering;
use std::collections::{HashMap, HashSet};
use std::io::{self, BufRead, Write};

const MAX_LOCI: usize = 6; // максимум генів в одній особині (3^6 = 729 станів)
const MAX_GENERATIONS: u32 = 50;

#[derive(Clone, Copy, PartialEq, Eq)]
enum Kind {
    Plant,
    Animal,
}

impl Kind {
    fn label(self) -> &'static str {
        match self {
            Kind::Plant => "рослини",
            Kind::Animal => "тварини",
        }
    }
    fn single(self) -> &'static str {
        match self {
            Kind::Plant => "рослину",
            Kind::Animal => "тварину",
        }
    }
}

#[derive(Default)]
struct Pool {
    loci: Vec<char>,                      // літери генів (великі), відсортовані
    traits: HashMap<char, (String, i32)>, // літера -> (назва, вага)
    inds: Vec<Vec<u8>>,                   // особини: для кожного гена 0=aa, 1=Aa, 2=AA
}

impl Pool {
    fn weight(&self, l: char) -> i32 {
        self.traits.get(&l).map(|t| t.1).unwrap_or(1)
    }

    fn trait_label(&self, l: char) -> String {
        match self.traits.get(&l) {
            Some((n, _)) if !n.is_empty() => format!("{} «{}»", l, n),
            _ => l.to_string(),
        }
    }

    fn add(&mut self, loci: Vec<char>, g: Vec<u8>) -> Result<(), String> {
        if self.inds.is_empty() {
            if loci.len() > MAX_LOCI {
                return Err(format!(
                    "Забагато генів ({}). Максимум — {}.",
                    loci.len(),
                    MAX_LOCI
                ));
            }
            self.loci = loci;
        } else if self.loci != loci {
            let s: String = self.loci.iter().collect();
            return Err(format!(
                "Набір генів має збігатися з уже введеними особинами: {}",
                s
            ));
        }
        self.inds.push(g);
        Ok(())
    }

    fn reset(&mut self) {
        self.inds.clear();
        self.loci.clear();
    }
}

#[derive(Default)]
struct Game {
    plants: Pool,
    animals: Pool,
}

impl Game {
    fn pool(&mut self, k: Kind) -> &mut Pool {
        match k {
            Kind::Plant => &mut self.plants,
            Kind::Animal => &mut self.animals,
        }
    }
}

// ───────────────────────── Робота з генотипами ─────────────────────────

fn lower(c: char) -> char {
    c.to_lowercase().next().unwrap_or(c)
}

fn upper(c: char) -> char {
    c.to_uppercase().next().unwrap_or(c)
}

/// "AaBbcc" -> (['A','B','C'], [1,1,0]); гени сортуються за літерою.
fn parse_genome(s: &str) -> Result<(Vec<char>, Vec<u8>), String> {
    let chars: Vec<char> = s
        .chars()
        .filter(|c| !c.is_whitespace() && *c != ',')
        .collect();
    if chars.is_empty() || chars.len() % 2 != 0 {
        return Err("Геном — це пари літер, наприклад: AaBbccDd".into());
    }
    let mut genes: Vec<(char, u8)> = Vec::new();
    for pair in chars.chunks(2) {
        let (a, b) = (pair[0], pair[1]);
        let cased = |c: char| c.is_uppercase() || c.is_lowercase();
        if !cased(a) || !cased(b) || lower(a) != lower(b) {
            return Err(format!(
                "Неправильна пара «{}{}»: обидві літери мають позначати один ген (AA, Aa або aa).",
                a, b
            ));
        }
        let letter = upper(a);
        if genes.iter().any(|g| g.0 == letter) {
            return Err(format!("Ген «{}» повторюється.", letter));
        }
        let dom = pair.iter().filter(|c| c.is_uppercase()).count() as u8;
        genes.push((letter, dom));
    }
    genes.sort_by_key(|g| g.0);
    Ok((
        genes.iter().map(|g| g.0).collect(),
        genes.iter().map(|g| g.1).collect(),
    ))
}

fn encode(g: &[u8]) -> usize {
    g.iter().rev().fold(0usize, |acc, &d| acc * 3 + d as usize)
}

fn digit(code: usize, i: usize) -> usize {
    code / 3usize.pow(i as u32) % 3
}

fn fmt_genotype(code: usize, loci: &[char]) -> String {
    let mut s = String::new();
    for (i, &l) in loci.iter().enumerate() {
        let lo = lower(l);
        match digit(code, i) {
            2 => {
                s.push(l);
                s.push(l);
            }
            1 => {
                s.push(l);
                s.push(lo);
            }
            _ => {
                s.push(lo);
                s.push(lo);
            }
        }
    }
    s
}

/// Рахунок особини: сума ваг генів, де проявляється домінантна ознака,
/// і "чистота" — скільки генів уже в бажаному гомозиготному стані.
fn evaluate(code: usize, pool: &Pool) -> (i32, usize) {
    let mut score = 0;
    let mut pure = 0;
    for (i, &l) in pool.loci.iter().enumerate() {
        let d = digit(code, i);
        let w = pool.weight(l);
        if d > 0 {
            score += w;
        }
        let want = if w >= 0 { 2 } else { 0 };
        if w == 0 || d == want {
            pure += 1;
        }
    }
    (score, pure)
}

fn phenotype(code: usize, pool: &Pool) -> String {
    let mut dom = Vec::new();
    let mut rec = Vec::new();
    for (i, &l) in pool.loci.iter().enumerate() {
        if digit(code, i) > 0 {
            dom.push(pool.trait_label(l));
        } else {
            rec.push(pool.trait_label(l));
        }
    }
    let j = |v: Vec<String>| if v.is_empty() { "—".to_string() } else { v.join(", ") };
    format!("домінантні: {}; рецесивні: {}", j(dom), j(rec))
}

/// Усі можливі нащадки схрещування a × b з імовірностями.
fn offspring_codes(a: usize, b: usize, n: usize) -> Vec<(usize, f64)> {
    let mut res = vec![(0usize, 1.0f64)];
    let (mut a, mut b) = (a, b);
    let mut mul = 1usize;
    for _ in 0..n {
        let x = a % 3;
        let y = b % 3;
        a /= 3;
        b /= 3;
        let pa = x as f64 / 2.0; // шанс передати домінантний алель
        let pb = y as f64 / 2.0;
        let d = [
            (1.0 - pa) * (1.0 - pb),
            pa * (1.0 - pb) + (1.0 - pa) * pb,
            pa * pb,
        ];
        let mut next = Vec::with_capacity(res.len() * 3);
        for &(c, p) in &res {
            for (k, &dk) in d.iter().enumerate() {
                if dk > 0.0 {
                    next.push((c + k * mul, p * dk));
                }
            }
        }
        res = next;
        mul *= 3;
    }
    res
}

// ───────────────────────── Пошук найкращого результату ─────────────────────────

#[derive(Clone)]
struct Recipe {
    level: u32,                    // у якому поколінні вперше досяжний
    parents: Option<(usize, usize)>, // None — це вихідна особина
    prob: f64,                     // шанс отримати цей генотип у такому схрещуванні
}

fn fmt_chance(p: f64) -> String {
    if p >= 0.9999 {
        return "100% — гарантовано".to_string();
    }
    format!("{:.2}% (≈ 1 з {} нащадків)", p * 100.0, (1.0 / p).round() as u64)
}

fn solve(pool: &Pool, kind: Kind) {
    if pool.inds.is_empty() {
        println!("Для {} ще немає особин. Додайте їх командою plant/animal.", kind.label());
        return;
    }
    let n = pool.loci.len();
    let total = 3usize.pow(n as u32);
    let mut rec: Vec<Option<Recipe>> = vec![None; total];
    let mut counts: HashMap<usize, u32> = HashMap::new();
    let mut known: Vec<usize> = Vec::new();

    for g in &pool.inds {
        let c = encode(g);
        *counts.entry(c).or_insert(0) += 1;
        if rec[c].is_none() {
            rec[c] = Some(Recipe { level: 0, parents: None, prob: 1.0 });
            known.push(c);
        }
    }

    // Пошук у ширину по поколіннях: у кожному поколінні схрещуємо нові генотипи з усіма відомими.
    let mut frontier = known.clone();
    let mut level = 0u32;
    while !frontier.is_empty() && level < MAX_GENERATIONS {
        level += 1;
        let mut in_front = vec![false; total];
        for &c in &frontier {
            in_front[c] = true;
        }
        let mut new_found: Vec<usize> = Vec::new();
        for &x in &frontier {
            for &y in &known {
                if in_front[y] && y < x {
                    continue; // пару (y, x) уже розглянули
                }
                if x == y {
                    // Рослина може запилювати сама себе; для тварин потрібні двоє однакових.
                    let ok = kind == Kind::Plant
                        || rec[x].as_ref().map_or(false, |r| r.level > 0)
                        || counts.get(&x).copied().unwrap_or(0) >= 2;
                    if !ok {
                        continue;
                    }
                }
                for (child, p) in offspring_codes(x, y, n) {
                    let action = match &rec[child] {
                        None => 1,
                        Some(r) if r.level == level && p > r.prob => 2,
                        _ => 0,
                    };
                    if action == 1 {
                        rec[child] = Some(Recipe { level, parents: Some((x, y)), prob: p });
                        new_found.push(child);
                    } else if action == 2 {
                        if let Some(r) = rec[child].as_mut() {
                            r.parents = Some((x, y));
                            r.prob = p;
                        }
                    }
                }
            }
        }
        known.extend(new_found.iter().copied());
        frontier = new_found;
    }

    let rc = |c: usize| rec[c].as_ref().unwrap();

    // Вибір найкращого: рахунок → чистота → менше поколінь → більший шанс.
    let mut best = known[0];
    for &c in &known[1..] {
        let (sc, pc) = evaluate(c, pool);
        let (sb, pb) = evaluate(best, pool);
        let ord = sc
            .cmp(&sb)
            .then(pc.cmp(&pb))
            .then(rc(best).level.cmp(&rc(c).level))
            .then(rc(c).prob.partial_cmp(&rc(best).prob).unwrap_or(Ordering::Equal));
        if ord == Ordering::Greater {
            best = c;
        }
    }

    let (score, pure) = evaluate(best, pool);
    let max_score: i32 = pool.loci.iter().map(|&l| pool.weight(l).max(0)).sum();

    println!("\n══════ Результат для {} ══════", kind.label());
    println!("Досяжних генотипів: {}.", known.len());
    println!("Найкращий досяжний генотип:  {}", fmt_genotype(best, &pool.loci));
    println!("  Ознаки:  {}", phenotype(best, pool));
    println!("  Рахунок: {} з можливих {}", score, max_score);
    if pure < n {
        println!(
            "  ⚠ Генотип не повністю чистий: його потомство розщеплюватиметься.\n    Це найкраще, що можна отримати з ваших особин."
        );
    } else {
        println!("  ✔ Чиста лінія: усі гени закріплені, потомство буде таким самим.");
    }

    if rc(best).level == 0 {
        println!("\nНайкраща особина вже є серед ваших — схрещування нічого не покращить.");
    } else {
        // Збираємо всі кроки, потрібні для найкращого генотипу.
        let mut needed: HashSet<usize> = HashSet::new();
        let mut stack = vec![best];
        while let Some(c) = stack.pop() {
            if !needed.insert(c) {
                continue;
            }
            if let Some((x, y)) = rc(c).parents {
                stack.push(x);
                stack.push(y);
            }
        }
        let mut originals: Vec<usize> =
            needed.iter().copied().filter(|&c| rc(c).level == 0).collect();
        originals.sort();
        let mut steps: Vec<usize> =
            needed.iter().copied().filter(|&c| rc(c).level > 0).collect();
        steps.sort_by_key(|&c| (rc(c).level, c));

        println!("\nВихідні особини, які потрібні:");
        for &c in &originals {
            println!("  • {}", fmt_genotype(c, &pool.loci));
        }

        println!("\nПлан схрещувань (поколінь: {}):", rc(best).level);
        let name = |c: usize| {
            let g = fmt_genotype(c, &pool.loci);
            if rc(c).level == 0 { format!("{} (ваш)", g) } else { g }
        };
        for (i, &c) in steps.iter().enumerate() {
            let r = rc(c);
            let (x, y) = r.parents.unwrap();
            let how = if x == y {
                match kind {
                    Kind::Plant => " [самозапилення]",
                    Kind::Animal => " [двоє однакових]",
                }
            } else {
                ""
            };
            println!(
                "  {:>2}. Покоління {}: {} × {}{}\n        → {}   шанс: {}",
                i + 1,
                r.level,
                name(x),
                name(y),
                how,
                fmt_genotype(c, &pool.loci),
                fmt_chance(r.prob)
            );
        }
        println!("\nПісля кожного кроку залишайте лише нащадків із потрібним генотипом.");
        if kind == Kind::Animal {
            println!("Для кожного схрещування потрібні самець і самка.");
        }
    }

    // Таблиця прогресу
    let max_level = known.iter().map(|&c| rc(c).level).max().unwrap_or(0);
    if max_level > 0 {
        println!("\nПрогрес по поколіннях (найкращий досяжний рахунок):");
        for lv in 0..=max_level {
            let mut top: Option<(usize, (i32, usize))> = None;
            for &c in &known {
                if rc(c).level > lv {
                    continue;
                }
                let e = evaluate(c, pool);
                if top.map_or(true, |(_, t)| e > t) {
                    top = Some((c, e));
                }
            }
            if let Some((c, e)) = top {
                println!(
                    "  Покоління {}: {}  рахунок {}",
                    lv,
                    fmt_genotype(c, &pool.loci),
                    e.0
                );
            }
        }
    }
    println!();
}

// ───────────────────────── Інтерфейс ─────────────────────────

fn help() {
    println!(
        "
Команди:
  plant <геном>     додати рослину, напр.  plant AaBbccDd      (скорочено: p)
  animal <геном>    додати тварину, напр.  animal MmWWss       (скорочено: a)
  trait <p|a> <літера> <вага> [назва]
                    цінність гена, напр.   trait p A 5 Врожайність
  list              показати введені особини
  del <p|a> <№>     видалити особину за номером
  clear <p|a|all>   очистити список
  solve [p|a]       знайти найкращий результат і план схрещувань
  demo              завантажити приклад
  help              ця довідка
  quit              вихід

Геном — пари літер: велика = домінантний алель, мала = рецесивний (AA, Aa, aa).
Порядок генів неважливий, але набір літер у всіх особин одного виду має збігатися.
Вага гена додається до рахунку, коли проявляється домінантна ознака
(від'ємна вага означає, що краще рецесивна). Вага за замовчуванням — 1.
"
    );
}

fn parse_kind(s: &str) -> Option<Kind> {
    match s.to_lowercase().as_str() {
        "p" | "plant" | "р" | "рослина" => Some(Kind::Plant),
        "a" | "animal" | "т" | "тварина" => Some(Kind::Animal),
        _ => None,
    }
}

fn add_genome(pool: &mut Pool, s: &str) -> Result<String, String> {
    let (loci, g) = parse_genome(s)?;
    let code = encode(&g);
    pool.add(loci, g)?;
    Ok(fmt_genotype(code, &pool.loci))
}

fn list_pool(pool: &Pool, kind: Kind) {
    println!("— {} —", kind.label().to_uppercase());
    if pool.inds.is_empty() {
        println!("  (порожньо)");
        return;
    }
    let genes: Vec<String> = pool
        .loci
        .iter()
        .map(|&l| format!("{} вага {}", pool.trait_label(l), pool.weight(l)))
        .collect();
    println!("  Гени: {}", genes.join("; "));
    for (i, g) in pool.inds.iter().enumerate() {
        let c = encode(g);
        println!(
            "  {}. {}   рахунок {}",
            i + 1,
            fmt_genotype(c, &pool.loci),
            evaluate(c, pool).0
        );
    }
}

fn load_demo(game: &mut Game) {
    game.plants = Pool::default();
    game.animals = Pool::default();
    for (l, n, w) in [
        ('A', "Врожайність", 5),
        ('B', "Стійкість до хвороб", 3),
        ('C', "Швидкість росту", 2),
        ('D', "Смак", 1),
    ] {
        game.plants.traits.insert(l, (n.to_string(), w));
    }
    for g in ["AaBbccDD", "aaBBCcdd", "AAbbCCDd"] {
        add_genome(&mut game.plants, g).unwrap();
    }
    for (l, n, w) in [('M', "Молочність", 4), ('W', "Якість вовни", 2), ('S', "Витривалість", 3)] {
        game.animals.traits.insert(l, (n.to_string(), w));
    }
    for g in ["MmWwss", "mmWWSs", "MMwwSS"] {
        add_genome(&mut game.animals, g).unwrap();
    }
    println!("Приклад завантажено. Спробуйте: list, потім solve.");
}

fn main() {
    println!("🌱🐄 Селекція рослин і тварин — введіть `help` для довідки, `demo` для прикладу.");
    let mut game = Game::default();
    let stdin = io::stdin();

    loop {
        print!("> ");
        io::stdout().flush().ok();
        let mut line = String::new();
        if stdin.lock().read_line(&mut line).unwrap_or(0) == 0 {
            break;
        }
        let line = line.trim();
        if line.is_empty() {
            continue;
        }
        let mut parts = line.split_whitespace();
        let cmd = parts.next().unwrap().to_lowercase();
        let rest: Vec<&str> = parts.collect();

        match cmd.as_str() {
            "help" | "?" => help(),
            "quit" | "exit" | "q" => break,
            "demo" => load_demo(&mut game),
            "plant" | "p" | "рослина" | "animal" | "a" | "тварина" => {
                let kind = parse_kind(&cmd).unwrap();
                match add_genome(game.pool(kind), &rest.join("")) {
                    Ok(g) => {
                        let n = game.pool(kind).inds.len();
                        println!("Додано {} №{}: {}", kind.single(), n, g);
                    }
                    Err(e) => println!("Помилка: {}", e),
                }
            }
            "trait" | "ген" => {
                if rest.len() < 3 {
                    println!("Формат: trait <p|a> <літера> <вага> [назва]");
                    continue;
                }
                let Some(kind) = parse_kind(rest[0]) else {
                    println!("Вкажіть p (рослини) або a (тварини).");
                    continue;
                };
                let mut chars = rest[1].chars();
                let (Some(l), None) = (chars.next(), chars.next()) else {
                    println!("Ген позначається однією літерою.");
                    continue;
                };
                let Ok(w) = rest[2].parse::<i32>() else {
                    println!("Вага має бути цілим числом.");
                    continue;
                };
                let l = upper(l);
                let pool = game.pool(kind);
                let prev = pool.traits.get(&l).map(|t| t.0.clone()).unwrap_or_default();
                let name = if rest.len() > 3 { rest[3..].join(" ") } else { prev };
                pool.traits.insert(l, (name, w));
                println!("Ген {} для {}: вага {}.", l, kind.label(), w);
            }
            "list" | "ls" => {
                list_pool(&game.plants, Kind::Plant);
                list_pool(&game.animals, Kind::Animal);
            }
            "del" => {
                if rest.len() < 2 {
                    println!("Формат: del <p|a> <№>");
                    continue;
                }
                let (Some(kind), Ok(i)) = (parse_kind(rest[0]), rest[1].parse::<usize>()) else {
                    println!("Формат: del <p|a> <№>");
                    continue;
                };
                let pool = game.pool(kind);
                if i == 0 || i > pool.inds.len() {
                    println!("Немає особини №{}.", i);
                } else {
                    pool.inds.remove(i - 1);
                    if pool.inds.is_empty() {
                        pool.reset();
                    }
                    println!("Видалено.");
                }
            }
            "clear" => match rest.first().map(|s| s.to_lowercase()) {
                Some(ref s) if s == "all" => {
                    game.plants.reset();
                    game.animals.reset();
                    println!("Усе очищено.");
                }
                Some(ref s) => match parse_kind(s) {
                    Some(k) => {
                        game.pool(k).reset();
                        println!("Очищено.");
                    }
                    None => println!("Формат: clear <p|a|all>"),
                },
                None => println!("Формат: clear <p|a|all>"),
            },
            "solve" | "s" => match rest.first() {
                Some(s) => match parse_kind(s) {
                    Some(Kind::Plant) => solve(&game.plants, Kind::Plant),
                    Some(Kind::Animal) => solve(&game.animals, Kind::Animal),
                    None => println!("Формат: solve [p|a]"),
                },
                None => {
                    if game.plants.inds.is_empty() && game.animals.inds.is_empty() {
                        println!("Спершу додайте особин (plant / animal) або введіть demo.");
                    }
                    if !game.plants.inds.is_empty() {
                        solve(&game.plants, Kind::Plant);
                    }
                    if !game.animals.inds.is_empty() {
                        solve(&game.animals, Kind::Animal);
                    }
                }
            },
            _ => println!("Невідома команда. Введіть help."),
        }
    }
}
