package main

import (
	"database/sql"
	"encoding/json"
	"errors"
	"fmt"
	"log"
	"net/http"
	"os"
	"strings"
	"time"
	"unicode"

	_ "modernc.org/sqlite"
)

type AnalyzeRequest struct {
	Sequence string `json:"sequence"`
}

type Analysis struct {
	ID                int64     `json:"id"`
	Sequence          string    `json:"sequence"`
	Length            int       `json:"length"`
	Counts            Counts    `json:"counts"`
	GCPercent         float64   `json:"gc_percent"`
	ReverseComplement string    `json:"reverse_complement"`
	MRNA              string    `json:"mrna"`
	Protein           string    `json:"protein"`
	ORFStart          int       `json:"orf_start,omitempty"`
	ORFEnd            int       `json:"orf_end,omitempty"`
	CreatedAt         time.Time `json:"created_at,omitempty"`
}

type Counts struct {
	A int `json:"A"`
	T int `json:"T"`
	G int `json:"G"`
	C int `json:"C"`
}

var codons = map[string]string{
	"UUU": "F", "UUC": "F", "UUA": "L", "UUG": "L", "UCU": "S", "UCC": "S", "UCA": "S", "UCG": "S",
	"UAU": "Y", "UAC": "Y", "UAA": "*", "UAG": "*", "UGU": "C", "UGC": "C", "UGA": "*", "UGG": "W",
	"CUU": "L", "CUC": "L", "CUA": "L", "CUG": "L", "CCU": "P", "CCC": "P", "CCA": "P", "CCG": "P",
	"CAU": "H", "CAC": "H", "CAA": "Q", "CAG": "Q", "CGU": "R", "CGC": "R", "CGA": "R", "CGG": "R",
	"AUU": "I", "AUC": "I", "AUA": "I", "AUG": "M", "ACU": "T", "ACC": "T", "ACA": "T", "ACG": "T",
	"AAU": "N", "AAC": "N", "AAA": "K", "AAG": "K", "AGU": "S", "AGC": "S", "AGA": "R", "AGG": "R",
	"GUU": "V", "GUC": "V", "GUA": "V", "GUG": "V", "GCU": "A", "GCC": "A", "GCA": "A", "GCG": "A",
	"GAU": "D", "GAC": "D", "GAA": "E", "GAG": "E", "GGU": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

func normalizeAndValidate(input string) (string, error) {
	var b strings.Builder
	for _, r := range input {
		if unicode.IsSpace(r) {
			continue
		}
		if r == '>' {
			return "", errors.New("FASTA-заголовки не поддерживаются: передайте только последовательность")
		}
		ch := unicode.ToUpper(r)
		if ch != 'A' && ch != 'T' && ch != 'G' && ch != 'C' {
			return "", fmt.Errorf("недопустимый символ %q: разрешены только A, T, G, C", r)
		}
		b.WriteRune(ch)
	}
	seq := b.String()
	if len(seq) == 0 {
		return "", errors.New("последовательность не должна быть пустой")
	}
	if len(seq) > 100000 {
		return "", errors.New("максимальная длина последовательности — 100000 нуклеотидов")
	}
	return seq, nil
}

func reverseComplement(seq string) string {
	comp := map[byte]byte{'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
	out := make([]byte, len(seq))
	for i := 0; i < len(seq); i++ {
		out[len(seq)-1-i] = comp[seq[i]]
	}
	return string(out)
}

func transcribe(seq string) string { return strings.ReplaceAll(seq, "T", "U") }

func translate(mrna string) string {
	var protein strings.Builder
	for i := 0; i+3 <= len(mrna); i += 3 {
		aa, ok := codons[mrna[i:i+3]]
		if !ok {
			break
		}
		if aa == "*" {
			break
		}
		protein.WriteString(aa)
	}
	return protein.String()
}

func analyze(seq string) Analysis {
	result := Analysis{Sequence: seq, Length: len(seq)}
	for i := range seq {
		switch seq[i] {
		case 'A':
			result.Counts.A++
		case 'T':
			result.Counts.T++
		case 'G':
			result.Counts.G++
		case 'C':
			result.Counts.C++
		}
	}
	if result.Length > 0 {
		result.GCPercent = float64(result.Counts.G+result.Counts.C) * 100 / float64(result.Length)
	}
	result.ReverseComplement = reverseComplement(seq)
	result.MRNA = transcribe(seq)
	result.Protein = translate(result.MRNA)

	start := strings.Index(seq, "ATG")
	if start >= 0 {
		for i := start; i+3 <= len(seq); i += 3 {
			codon := seq[i : i+3]
			if i > start && (codon == "TAA" || codon == "TAG" || codon == "TGA") {
				result.ORFStart = start + 1
				result.ORFEnd = i + 3
				break
			}
		}
	}
	return result
}

func main() {
	dbPath := os.Getenv("BIOLAB_DB")
	if dbPath == "" {
		dbPath = "biolab.db"
	}
	db, err := sql.Open("sqlite", dbPath)
	if err != nil {
		log.Fatal(err)
	}
	defer db.Close()
	db.SetMaxOpenConns(1)
	if _, err := db.Exec(`PRAGMA journal_mode=WAL; CREATE TABLE IF NOT EXISTS analyses (
		id INTEGER PRIMARY KEY AUTOINCREMENT,
		sequence TEXT NOT NULL,
		length INTEGER NOT NULL,
		counts_json TEXT NOT NULL,
		gc_percent REAL NOT NULL,
		reverse_complement TEXT NOT NULL,
		mrna TEXT NOT NULL,
		protein TEXT NOT NULL,
		orf_start INTEGER NOT NULL,
		orf_end INTEGER NOT NULL,
		created_at TEXT NOT NULL
	)`); err != nil {
		log.Fatal(err)
	}

	mux := http.NewServeMux()
	mux.HandleFunc("/api/health", func(w http.ResponseWriter, r *http.Request) {
		if r.Method != http.MethodGet {
			methodNotAllowed(w)
			return
		}
		writeJSON(w, http.StatusOK, map[string]string{"status": "ok"})
	})
	mux.HandleFunc("/api/analyze", func(w http.ResponseWriter, r *http.Request) {
		if r.Method != http.MethodPost {
			methodNotAllowed(w)
			return
		}
		r.Body = http.MaxBytesReader(w, r.Body, 1<<20)
		var req AnalyzeRequest
		dec := json.NewDecoder(r.Body)
		dec.DisallowUnknownFields()
		if err := dec.Decode(&req); err != nil {
			writeJSON(w, http.StatusBadRequest, map[string]string{"error": "Ожидается JSON с полем sequence"})
			return
		}
		seq, err := normalizeAndValidate(req.Sequence)
		if err != nil {
			writeJSON(w, http.StatusBadRequest, map[string]string{"error": err.Error()})
			return
		}
		result := analyze(seq)
		countsJSON, _ := json.Marshal(result.Counts)
		now := time.Now().UTC()
		res, err := db.Exec(`INSERT INTO analyses(sequence,length,counts_json,gc_percent,reverse_complement,mrna,protein,orf_start,orf_end,created_at)
			VALUES(?,?,?,?,?,?,?,?,?,?)`, result.Sequence, result.Length, string(countsJSON), result.GCPercent, result.ReverseComplement, result.MRNA, result.Protein, result.ORFStart, result.ORFEnd, now.Format(time.RFC3339))
		if err != nil {
			writeJSON(w, http.StatusInternalServerError, map[string]string{"error": "Не удалось сохранить анализ"})
			return
		}
		result.ID, _ = res.LastInsertId()
		result.CreatedAt = now
		writeJSON(w, http.StatusOK, result)
	})
	mux.HandleFunc("/api/history", func(w http.ResponseWriter, r *http.Request) {
		if r.Method != http.MethodGet {
			methodNotAllowed(w)
			return
		}
		rows, err := db.Query(`SELECT id,sequence,length,counts_json,gc_percent,reverse_complement,mrna,protein,orf_start,orf_end,created_at FROM analyses ORDER BY id DESC LIMIT 50`)
		if err != nil {
			writeJSON(w, http.StatusInternalServerError, map[string]string{"error": "Не удалось прочитать историю"})
			return
		}
		defer rows.Close()
		history := make([]Analysis, 0)
		for rows.Next() {
			var a Analysis
			var countsJSON, created string
			if err := rows.Scan(&a.ID, &a.Sequence, &a.Length, &countsJSON, &a.GCPercent, &a.ReverseComplement, &a.MRNA, &a.Protein, &a.ORFStart, &a.ORFEnd, &created); err != nil {
				writeJSON(w, http.StatusInternalServerError, map[string]string{"error": "Ошибка чтения истории"})
				return
			}
			_ = json.Unmarshal([]byte(countsJSON), &a.Counts)
			a.CreatedAt, _ = time.Parse(time.RFC3339, created)
			history = append(history, a)
		}
		writeJSON(w, http.StatusOK, history)
	})
	mux.Handle("/", http.FileServer(http.Dir("static")))

	addr := ":8080"
	log.Printf("BioLab запущен: http://localhost%s", addr)
	log.Fatal(http.ListenAndServe(addr, withCORS(mux)))
}

func methodNotAllowed(w http.ResponseWriter) {
	writeJSON(w, http.StatusMethodNotAllowed, map[string]string{"error": "Метод не поддерживается"})
}
func writeJSON(w http.ResponseWriter, status int, value any) {
	w.Header().Set("Content-Type", "application/json; charset=utf-8")
	w.WriteHeader(status)
	_ = json.NewEncoder(w).Encode(value)
}
func withCORS(next http.Handler) http.Handler {
	return http.HandlerFunc(func(w http.ResponseWriter, r *http.Request) {
		w.Header().Set("Access-Control-Allow-Origin", "*")
		w.Header().Set("Access-Control-Allow-Headers", "Content-Type")
		w.Header().Set("Access-Control-Allow-Methods", "GET, POST, OPTIONS")
		if r.Method == http.MethodOptions {
			w.WriteHeader(http.StatusNoContent)
			return
		}
		next.ServeHTTP(w, r)
	})
}
