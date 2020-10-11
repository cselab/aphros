// Created by Petr Karnakov on 11.10.2020
// Copyright 2020 ETH Zurich

#include "argparse.h"

struct ArgumentParser::Imp {
  Imp(ArgumentParser& owner, std::string desc, bool isroot)
      : owner_(owner), desc_(desc), isroot_(isroot) {
    AddSwitch({"--help", "-h"}).Help("Print help and exit");
  }
  Proxy AddSwitch(const std::vector<std::string>& names) {
    fassert(!names.empty(), "Empty list of names");
    const std::string key = StripDash(names[0]);
    const bool optional = IsOptional(names[0]);
    fassert(
        optional, "Switch name must start with '-', got '" + names[0] + "'");
    CheckNameConsistency(names);
    Entry entry{key, optional, true, 0};
    for (auto name : names) {
      CheckNameConflict(name);
      entries_[name] = entry;
      names_[key].push_back(name);
    }
    entries_[key] = entry;
    known_args_.Int.Set(key, 0);
    opt_keys_.push_back(key);
    return Proxy(owner_, key);
  }

  template <class T>
  Proxy AddVariable(
      const std::vector<std::string>& names, T defaultval, bool hasdefault) {
    fassert(!names.empty(), "Empty list of names");
    const std::string key = StripDash(names[0]);
    const bool optional = IsOptional(names[0]);
    CheckNameConsistency(names);
    Entry entry{key, optional, hasdefault, 1};
    for (auto name : names) {
      CheckNameConflict(name);
      entries_[name] = entry;
      names_[key].push_back(name);
    }
    entries_[key] = entry;
    known_args_.template Get<T>().Set(key, defaultval);
    if (optional) {
      opt_keys_.push_back(key);
    } else if (!optional) {
      pos_keys_.push_back(key);
    }
    return Proxy(owner_, key);
  }

  Vars ParseArgs(
      std::vector<std::string> argv, std::string program) const {
    Vars args;

    args.Int.Set("FAIL", 0);

    auto rassert = [&](bool cond, std::string msg) {
      if (!cond) {
        if (isroot_) {
          std::cerr << msg << '\n';
        }
        args.Int["FAIL"] = 1;
      }
      return cond;
    };

    for (auto& it : entries_) {
      auto entry = it.second;
      if (entry.hasdefault) {
        auto key = entry.key;
        auto type = known_args_.GetTypeName(key);
        args.SetStr(type, key, known_args_.GetStr(type, key));
      }
    }

    enum class S { name, optional, positional };
    S s = S::name;
    Entry entry;
    size_t ipos = 0; // index in pos_keys_
    for (size_t i = 0; i < argv.size(); ++i) {
      auto str = argv[i];
      switch (s) {
        case S::name:
          if (IsOptional(str)) {
            if (rassert(entries_.count(str), "Unknown argument: " + str)) {
              entry = entries_.at(str);
              fassert(entry.nargs == 0 || entry.nargs == 1);
              if (entry.nargs == 1) {
                s = S::optional;
              } else { // expecting a switch, type int
                fassert(
                    known_args_.GetTypeName(entry.key) == "int",
                    "Expected type int for name='" + str + "'");
                args.SetStr(known_args_.GetTypeName(entry.key), entry.key, "1");
              }
            }
          } else { // positional
            s = S::positional;
            --i;
          }
          break;
        case S::optional:
          args.SetStr(known_args_.GetTypeName(entry.key), entry.key, str);
          s = S::name;
          break;
        case S::positional:
          if (rassert(
                  ipos < pos_keys_.size(),
                  util::Format(
                      "Too many positional arguments: {}, expected {}",
                      ipos + 1, pos_keys_.size()))) {
            entry = entries_.at(pos_keys_[ipos]);
            args.SetStr(known_args_.GetTypeName(entry.key), entry.key, str);
          }
          ++ipos;
          s = S::name;
          break;
      }
    }

    if (args.Int["help"]) {
      if (isroot_) {
        PrintHelp(std::cerr, true, program);
      }
      args.Int.Set("EXIT", 0);
      return args;
    }

    for (; ipos < pos_keys_.size(); ++ipos) {
      rassert(
          entries_.at(pos_keys_[ipos]).hasdefault,
          "Missing value for positional argument " + ToUpper(pos_keys_[ipos]) +
              " without a default");
    }

    if (args.Int["FAIL"]) {
      if (isroot_) {
        PrintHelp(std::cerr, false, program);
      }
      args.Int.Set("EXIT", 1);
    }

    return args;
  }

  void PrintHelp(std::ostream& out, bool full, std::string program) const {
    out << "usage: " << program;
    for (auto key : opt_keys_) {
      out << " [";
      bool first = true;
      for (auto name : names_.at(key)) {
        out << (first ? "" : "|") << name;
        first = false;
      }
      if (entries_.at(key).nargs == 1) {
        out << ' ' << ToUpper(key);
      }
      out << "]";
    }
    for (auto key : pos_keys_) {
      out << " " << ToUpper(key);
    }
    out << std::endl;

    if (full) {
      if (desc_.length()) {
        out << '\n' << desc_ << '\n';
      }

      out << "\npositional arguments:\n";
      for (auto key : pos_keys_) {
        out << PadRight(ToUpper(key), 20);
        auto help = help_.count(key) ? help_.at(key) : "";
        out << help;
        auto entry = entries_.at(key);
        if (entry.hasdefault && entry.nargs) {
          out << (help.length() ? ". " : "") << "Default is "
              << known_args_.GetStr(known_args_.GetTypeName(key), key);
        }
        out << '\n';
      }

      out << "\noptional arguments:\n";
      for (auto key : opt_keys_) {
        std::string pre;
        bool first = true;
        for (auto name : names_.at(key)) {
          pre += (first ? "" : ", ") + name;
          first = false;
        }
        out << PadRight(pre, 20);
        auto help = help_.count(key) ? help_.at(key) : "";
        out << help;
        auto entry = entries_.at(key);
        if (entry.hasdefault && entry.nargs) {
          out << (help.length() ? ". " : "") << "Default is "
              << known_args_.GetStr(known_args_.GetTypeName(key), key);
        }
        out << '\n';
      }
    }
  }

  struct Entry {
    std::string key; // target key in vars
    bool optional; // if arguments starts with '-'
    bool hasdefault; // if argument has default value
    int nargs; // number of arguments to read,
               // values larger than 1 only for std::vector<double>
               // -1 to keep reading while numbers
  };
  static std::string StripDash(std::string name) {
    size_t i = 0;
    while (i < name.length() && name[i] == '-') {
      ++i;
    }
    return name.substr(i);
  }
  static bool IsNumber(std::string str) {
    std::stringstream buf(str);
    double a;
    buf >> a;
    return !buf.fail();
  }
  static bool StartsWithDash(std::string name) {
    return name.length() && name[0] == '-';
  }
  static bool IsOptional(std::string name) {
    return StartsWithDash(name) && !IsNumber(name);
  }
  static void CheckNameConsistency(std::vector<std::string> names) {
    if (names.empty()) {
      return;
    }
    for (auto& name : names) {
      fassert(
          IsOptional(name) == IsOptional(names[0]),
          "Names have different types (optional or positional): " + name +
              " and " + names[0]);
    }
  }
  void CheckNameConflict(std::string name) const {
    auto it = entries_.find(name);
    fassert(
        it == entries_.end(),
        "Name '" + name + "' already used with key '" + it->second.key + "'");
  }
  static std::string ToUpper(std::string s) {
    for (auto& c : s) {
      c = std::toupper(c);
    }
    return s;
  }
  static std::string PadRight(std::string s, int width) {
    return s + std::string(std::max<int>(0, width - s.length()), ' ');
  }
  
  ArgumentParser& owner_;
  std::string desc_;
  const bool isroot_;
  std::map<std::string, Entry> entries_; // name or key to entry
  Vars known_args_;
  std::vector<std::string> pos_keys_;
  std::vector<std::string> opt_keys_;
  std::map<std::string, std::string> help_; // key to help
  std::map<std::string, std::vector<std::string>> names_; // key to names
};

ArgumentParser::Proxy::Proxy(ArgumentParser& parser, std::string key)
    : parser_(parser), key_(key) {}

ArgumentParser::Proxy ArgumentParser::Proxy::Help(std::string help) {
  parser_.imp->help_[key_] = help;
  return *this;
}

ArgumentParser::ArgumentParser(std::string desc, bool isroot)
    : imp(new Imp(*this, desc, isroot)) {}

ArgumentParser::~ArgumentParser() = default;

auto ArgumentParser::AddSwitch(const std::vector<std::string>& names) -> Proxy {
  return imp->AddSwitch(names);
}

template <class T>
auto ArgumentParser::AddVariable(
    const std::vector<std::string>& names, T defaultval, bool hasdefault)
    -> Proxy {
  return imp->AddVariable(names, defaultval, hasdefault);
}

template auto ArgumentParser::AddVariable(
    const std::vector<std::string>& names, double defaultval,
    bool hasdefault)->Proxy;
template auto ArgumentParser::AddVariable(
    const std::vector<std::string>& names, int defaultval,
    bool hasdefault)->Proxy;
template auto ArgumentParser::AddVariable(
    const std::vector<std::string>& names, std::string defaultval,
    bool hasdefault)->Proxy;
template auto ArgumentParser::AddVariable(
    const std::vector<std::string>& names, std::vector<double> defaultval,
    bool hasdefault)->Proxy;

auto ArgumentParser::ParseArgs(
    std::vector<std::string> argv, std::string program) const -> Vars {
  return imp->ParseArgs(argv, program);
}

auto ArgumentParser::GetKnownArgs() const -> const Vars& {
  return imp->known_args_;
}

void ArgumentParser::PrintHelp(
    std::ostream& out, bool full, std::string program) const {
  imp->PrintHelp(out, full, program);
}
