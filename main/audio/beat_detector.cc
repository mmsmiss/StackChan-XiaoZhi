#include "beat_detector.h"

#include <cmath>
#include <algorithm>
#include <esp_timer.h>
#include <esp_log.h>

static const char* TAG = "BeatDetector";

namespace {

constexpr float kPi = 3.14159265358979323846f;

// 低频鼓点检测频段：bin 1 ~ 8（16kHz / 512 点 = 31.25Hz/格，约 31~250Hz，覆盖鼓点/贝斯）
constexpr int kLowBinStart = 1;
constexpr int kLowBinEnd = 8;

// 节拍判定阈值
constexpr float kBeatAbsMinEnergy = 0.05f;   // 低频绝对能量下限（归一化幅度）
constexpr float kBeatRatio = 1.5f;            // 低频能量相对历史均值突增倍数
constexpr int kBeatMinGapMs = 240;            // 最小节拍间隔（约 250 BPM 上限）

// 迭代式 radix-2 复数 FFT（就地，长度须为 2 的幂）
void Fft(std::vector<float>& re, std::vector<float>& im) {
    const size_t n = re.size();
    // 位反转
    for (size_t i = 1, j = 0; i < n; ++i) {
        size_t bit = n >> 1;
        for (; j & bit; bit >>= 1) j ^= bit;
        j ^= bit;
        if (i < j) {
            std::swap(re[i], re[j]);
            std::swap(im[i], im[j]);
        }
    }
    // 蝶形
    for (size_t len = 2; len <= n; len <<= 1) {
        float ang = -2.0f * kPi / len;
        float w_re = std::cos(ang), w_im = std::sin(ang);
        for (size_t i = 0; i < n; i += len) {
            float cur_re = 1.0f, cur_im = 0.0f;
            for (size_t j = 0; j < len / 2; ++j) {
                size_t a = i + j, b = i + j + len / 2;
                float t_re = re[b] * cur_re - im[b] * cur_im;
                float t_im = re[b] * cur_im + im[b] * cur_re;
                re[b] = re[a] - t_re;
                im[b] = im[a] - t_im;
                re[a] += t_re;
                im[a] += t_im;
                float nw_re = cur_re * w_re - cur_im * w_im;
                cur_im = cur_re * w_im + cur_im * w_re;
                cur_re = nw_re;
            }
        }
    }
}

}  // namespace

BeatDetector::BeatDetector() {
    buffer_.resize(kFftSize, 0.0f);
    // Hann 窗
    window_.resize(kFftSize);
    for (uint16_t i = 0; i < kFftSize; ++i) {
        window_[i] = 0.5f - 0.5f * std::cos(2.0f * kPi * i / (kFftSize - 1));
    }
}

void BeatDetector::Reset() {
    buffer_pos_ = 0;
    last_beat_ms_ = 0;
    avg_low_energy_ = 0.0f;
    std::fill(buffer_.begin(), buffer_.end(), 0.0f);
}

bool BeatDetector::Process(const std::vector<int16_t>& pcm) {
    if (pcm.empty()) return false;

    for (int16_t s : pcm) {
        buffer_[buffer_pos_++] = static_cast<float>(s) / 32768.0f;
        if (buffer_pos_ >= kFftSize) {
            buffer_pos_ = 0;
            int64_t now_ms = esp_timer_get_time() / 1000;
            if (ComputeBeat(buffer_, now_ms)) {
                return true;
            }
        }
    }
    return false;
}

bool BeatDetector::ComputeBeat(const std::vector<float>& samples, int64_t now_ms) {
    // 最小间隔限制
    if (now_ms - last_beat_ms_ < kBeatMinGapMs) {
        return false;
    }

    std::vector<float> re(kFftSize), im(kFftSize, 0.0f);
    for (uint16_t i = 0; i < kFftSize; ++i) {
        re[i] = samples[i] * window_[i];
    }
    Fft(re, im);

    // 低频段能量（幅度谱，只用前半）
    float low_energy = 0.0f;
    for (int b = kLowBinStart; b <= kLowBinEnd && b < kFftSize / 2; ++b) {
        low_energy += std::sqrt(re[b] * re[b] + im[b] * im[b]);
    }
    low_energy /= (kLowBinEnd - kLowBinStart + 1);  // 每格平均

    // 滑动均值（慢速，跟随背景音量）
    avg_low_energy_ = avg_low_energy_ * 0.95f + low_energy * 0.05f;

    // 节拍判定：低频能量明显高于历史均值（突增 = 鼓点）且超过绝对下限
    bool is_beat = (low_energy > kBeatAbsMinEnergy * kFftSize * 0.5f) &&
                   (low_energy > avg_low_energy_ * kBeatRatio);
    if (is_beat) {
        last_beat_ms_ = now_ms;
        return true;
    }
    return false;
}
