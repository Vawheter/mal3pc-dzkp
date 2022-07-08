# mal3pc-dzkp

### 安装Rust
参考[官方](https://www.rust-lang.org/zh-CN/tools/install)推荐安装方式
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh

### 运行测试
测试函数一般都在mod.rs

##### 运行某个测试
切换到dzkp目录，运行

cargo test xxx --release -- --nocapture

xxx是测试函数名称，--release开启优化，-- --nocapture打印输出

##### 运行所有测试
cargo test --release -- --nocapture