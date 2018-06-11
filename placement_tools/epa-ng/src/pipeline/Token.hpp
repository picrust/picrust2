#pragma once

#include <cereal/types/base_class.hpp>

enum class token_status {DATA, END};

class Token
{
public:

  Token()   = default;
  ~Token()  = default;

  virtual bool valid() const final
  {
    return (status_ != token_status::END);
  }

  virtual void status(const token_status& s) final
  {
    status_ = s;
  }

  virtual token_status status() const final
  {
    return status_;
  }

  virtual void is_last(const bool b) final
  {
    if (b) {
      status_ = token_status::END;
    }
  }

  template <class Archive>
  void serialize( Archive & ar )
  { ar( status_ ); }

private:

  token_status status_ = token_status::DATA;
  
};

class VoidToken : public Token
{
public:
  VoidToken() = default;
  ~VoidToken()= default;

  size_t size() {return 0;}

  void clear () {;}

  template <class Archive>
  void serialize( Archive & ar )
  { ar( *static_cast<Token*>( this ) ); }
};
